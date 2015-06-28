/* Copyright (c) 2015, The University of Oxford
   BDS 3-Clause license (see LICENSE) */
#include <sys/time.h>
#include <casacore/tables/Tables.h>
#include <casacore/ms/MeasurementSets.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <string>
#include <complex>
#include <cmath>
#include <vector>
#include <cassert>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::complex;
using std::flush;

using casa::MeasurementSet;
using casa::Table;
using casa::TableColumn;
using casa::TableExprNode;
using casa::Array;
using casa::Matrix;
using casa::ArrayColumn;
using casa::ScalarColumn;
using casa::ROMSColumns;
using casa::IPosition;
using casa::Double;
using casa::Float;
using casa::Complex;
using casa::Vector;
using casa::uInt;
using casa::Int;
using casa::Bool;

double timer_elapsed(const struct timeval& t1, const struct timeval& t2) {
    return (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1.0e6;
}

void copy_ms(string ms_in, string ms_out) {
    const Table ms(ms_in);
    ms.deepCopy(ms_out, Table::New, false, Table::LittleEndian, true);
    // TODO(BM) get list of sub-tables automatically?
    vector<string> sub_tables;
    sub_tables.push_back("/ANTENNA");
    sub_tables.push_back("/OBSERVATION");
    sub_tables.push_back("/SPECTRAL_WINDOW");
    sub_tables.push_back("/DATA_DESCRIPTION");
    sub_tables.push_back("/FEED");
    sub_tables.push_back("/FIELD");
    sub_tables.push_back("/HISTORY");
    sub_tables.push_back("/POLARIZATION");
    for (vector<string>::iterator it = sub_tables.begin();
                    it != sub_tables.end(); ++it) {
        const Table sub_table(ms_in+*it);
        sub_table.deepCopy(ms_out+*it, Table::New, true,
                           Table::LittleEndian, false);
    }
}


/* arcsinc(x) function taken from Obit. Uses Newton-Raphson method.  */
double inv_sinc(double value) {
    double x1 = 0.001;
    for (int i = 0; i < 1000; ++i) {
        double x0 = x1;
        double a = x0 * M_PI;
        x1 = x0 - ((sin(a) / a) - value) /
                        ((a * cos(a) - M_PI * sin(a)) / (a * a));
        if (fabs(x1 - x0) < 1.0e-6)
            break;
    }
    return x1;
}


/*
 * TODO(BM):
 * - Load and average in blocks (base block length on dt_max).
 * - Exploit extra performance of blocks and addrow(init=false)
 * - Handle polarisation
 * - Handle MODEL_DATA and CORRECTED_DATA columns.
 */
int main(int argc, char** argv) {
    if (argc != 5) {
        cout << "Usage:" << endl;
        cout << "  ./bda <ms> <max_fact> <fov> <idt_max>" << endl << endl;
        cout << "Where:" << endl;
        cout << "  ms       Path of measurement set to be averaged." << endl;
        cout << "  max_fact Maximum amplitude reduction factor. (eg. 1.01)\n";
        cout << "  fov      Field-of-view radius in degrees associated\n";
        cout << "           with the amplitude reduction factor" << endl;
        cout << "  idt_max  Maximum averaging time in terms of non-averaged\n";
        cout << "           correlator dumps." << endl << endl;
        cout << "Example:" << endl;
        cout << "  ./bda vis.ms 1.01 0.9 4" << endl;
        return 0;
    }
    string filename = string(argv[1]);
    double max_fact = atof(argv[2]);
    double fov = atof(argv[3]);
    int dt_max_i = atoi(argv[4]);
    string rootname = filename.substr(0, filename.find_last_of("."));
    string filename_out = rootname+"_bda.ms";

    // Create a total run time timer.
    struct timeval timer_all[2];
    gettimeofday(&timer_all[0], NULL);

    // Copy input measurement set.
    struct timeval timer_copy[2];
    gettimeofday(&timer_copy[0], NULL);
    cout << "+ Copying data ... " << flush;
    copy_ms(filename, filename_out);
    gettimeofday(&timer_copy[1], NULL);
    printf("\r+ Copying data completed in %.3f s\n",
           timer_elapsed(timer_copy[0], timer_copy[1]));

    // Load input measurement set.
    // ************************************************************************
    struct timeval timer_load[2];
    gettimeofday(&timer_load[0], NULL);
    cout << "+ Loading data ..." << flush;
    const Table ms(filename);
    uInt num_rows = ms.nrow();

    // Load antenna table to get number of antennas.
    cout << "\r" "+ Loading data : Antenna table." << flush;
    const Table ant(filename+"/ANTENNA");
    uInt num_antennas = ant.nrow();

    // Work out number of baselines and number of times.
    uInt num_baselines = num_antennas * (num_antennas - 1) / 2;
    uInt num_times = num_rows / num_baselines;

    // Load observation table to get time range (needed for delta_t).
    cout << "\r" "+ Loading data : Observation table." << flush;
    const Table obs(filename+"/OBSERVATION");
    const ArrayColumn<Double> time_range(obs, "TIME_RANGE");
    const Array<Double> time_range_cell = time_range.get(0);
    double t0 = time_range_cell.data()[0];
    double t1 = time_range_cell.data()[1];
    double obs_length_s = t1 - t0;
    double delta_t = obs_length_s / num_times;

    // Load frequency from MS. (require 1 freq)
    cout << "\r" "+ Loading data : Spectral window table." << flush;
    const Table spectral_window(filename+"/SPECTRAL_WINDOW");
    const ArrayColumn<Double> col_chan_freq(spectral_window, "CHAN_FREQ");
    const Vector<Double> freqs0 = col_chan_freq.get(0);
    int num_freqs = freqs0.size();
    if (num_freqs != 1) {
        cerr << endl << endl;
        cerr << "ERROR: This function only works for MS with one "
                        "frequency! Detected: " << num_freqs << endl;
        return 1;
    }
    double freq = freqs0.data()[0];
    double wavelength = 299792458.0 / freq;

    // Load observation table to get time range (needed for delta_t).
    cout << "\r" "+ Loading data : Polarisation table." << flush;
    const Table pol(filename+"/POLARIZATION");
    const ArrayColumn<Int> corr_type(pol, "CORR_TYPE");
    const Vector<Int> corr_type_cell = corr_type.get(0);
    int num_pols = corr_type_cell.size();

    cout << "\r" "+ Loading data : Main table columns.        " << flush;
    const ArrayColumn<Double> col_uvw(ms, "UVW");
    const ArrayColumn<Complex> col_data(ms, "DATA");
    const ScalarColumn<Double> col_time(ms, "TIME");
    const Array<Double> uvw_ = col_uvw.getColumn();
    const Array<Complex> data_ = col_data.getColumn();
    const Vector<Double> time_ = col_time.getColumn();
    const Complex* data = data_.data();
    const double* time = time_.data();
    const double* uvw = uvw_.data();
    vector<double> uu_(num_rows);
    vector<double> vv_(num_rows);
    vector<double> ww_(num_rows);
    double* uu = &uu_[0];
    double* vv = &vv_[0];
    double* ww = &ww_[0];
    // Make sure the data column dimension matches the number of polarisations.
    assert(col_data.shapeColumn()[0] == num_pols);
    assert(col_data.shapeColumn()[1] == 1);

    cout << "\r" "+ Loading data starting rearrange of uvw ... " << flush;
    for (uInt i = 0; i < num_rows; ++i) {
        uu[i] = uvw[i*3 + 0];
        vv[i] = uvw[i*3 + 1];
        ww[i] = uvw[i*3 + 2];
    }
    gettimeofday(&timer_load[1], NULL);
    printf("\r+ Loading data completed in %.3f s                 \n",
           timer_elapsed(timer_load[0], timer_load[1]));
    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------

    // Evaluate averaging parameters.
    double dt_max = dt_max_i * delta_t;
    double delta_uvw = inv_sinc(1.0 / max_fact) / (fov * (M_PI / 180.0));
    delta_uvw *= wavelength;
    double duvw_max = delta_uvw;

    printf("%s\n", string(80, '-').c_str());
    printf("Input MS        : %s\n", filename.c_str());
    printf("  No. rows      : %i\n", num_rows);
    printf("  No. antennas  : %i\n", num_antennas);
    printf("  No. baselines : %i\n", num_baselines);
    printf("  No. times     : %i\n", num_times);
    printf("  No. pols      : %i\n", num_pols);
    printf("  freq.         : %.16f\n", freq);
    printf("  delta_t       : %.16f\n", delta_t);
    printf("Averaging:\n");
    printf("  MS out        : %s\n", filename_out.c_str());
    printf("  max_fact      : %f\n", max_fact);
    printf("  FoV           : %f deg.\n", fov);
    printf("  dt_max_i      : %i\n", dt_max_i);
    printf("  dt_max        : %f s\n", dt_max);
    printf("  duvw_max      : %f m\n", duvw_max);
    printf("%s\n", string(80, '-').c_str());

    cout << "+ Opening output table..." << flush;
    // Open the output MS and get references to columns.
    Table ms_out(filename_out, Table::Update);
    ArrayColumn<Double> col_uvw_out(ms_out, "UVW");
    ScalarColumn<Double> col_time_out(ms_out, "TIME");
    ScalarColumn<Double> col_time_centroid_out(ms_out, "TIME_CENTROID");
    ScalarColumn<Int> col_ant1_out(ms_out, "ANTENNA1");
    ScalarColumn<Int> col_ant2_out(ms_out, "ANTENNA2");
    ArrayColumn<Bool> col_flag_out(ms_out, "FLAG");
    ScalarColumn<Double> col_interval_out(ms_out, "INTERVAL");
    ScalarColumn<Double> col_exposure_out(ms_out, "EXPOSURE");
    ArrayColumn<Complex> col_data_out(ms_out, "DATA");
    ArrayColumn<Float> col_weight_out(ms_out, "WEIGHT");
    ArrayColumn<Float> col_sigma_out(ms_out, "SIGMA");
    cout << "\r" "+ Output table and columns ready for writing." << endl;

    // Buffers of deltas along the baseline in the current average.
    std::vector<double> duvw_(num_baselines, 0.0);
    std::vector<double> dt_(num_baselines, 0);
    double* duvw = &duvw_[0];
    double* dt = &dt_[0];
    // Buffers of the current average along the baseline.
    vector<double> ave_uu_(num_baselines, 0.0);
    vector<double> ave_vv_(num_baselines, 0.0);
    vector<double> ave_ww_(num_baselines, 0.0);
    vector<double> ave_t_(num_baselines, 0.0);
    vector<Complex> ave_data_(num_baselines * num_pols, 0.0);

    vector<uInt> ave_count_(num_baselines, 0);
    double* ave_uu = &ave_uu_[0];
    double* ave_vv = &ave_vv_[0];
    double* ave_ww = &ave_ww_[0];
    double* ave_t = &ave_t_[0];
    uInt* ave_count = &ave_count_[0];

    Vector<Double> bda_uvw(3);
    Matrix<Complex> bda_data(num_pols, 1);
    Matrix<Bool> flag(num_pols, num_freqs, false);
    Complex* bda_data_p = bda_data[0].data();
    Complex* ave_data = &ave_data_[0];

    struct timeval timer_loop[2];
    gettimeofday(&timer_loop[0], 0);
    int bda_row = 0;
    for (uInt t = 0; t < num_times; ++t) {
        for (uInt b = 0, a1 = 0; a1 < num_antennas; ++a1) {
            for (uInt a2 = a1+1; a2 < num_antennas; ++a2, ++b) {
                uInt row = t * num_baselines + b;

                // Get coordinates and time for the current baseline.
                double b_uu = uu[row];
                double b_vv = vv[row];
                double b_ww = ww[row];
                double b_t  = time[row];

                // Accumulate into the average.
                ave_uu[b] += b_uu;
                ave_vv[b] += b_vv;
                ave_ww[b] += b_ww;
                ave_t[b] += b_t;
                for (int p = 0; p < num_pols; ++p) {
                    ave_data[b * num_pols + p] += data[row * num_pols + p];
                }

                ave_count[b] += 1;

                // Look ahead to the next next point on the baseline and see if
                // this is also in the average.
                double b_duvw = 0.0;
                double b_dt = 0.0;
                if (t < num_times - 1) {
                    uInt row_next = (t + 1) * num_baselines + b;
                    double b_uu1 = uu[row_next];
                    double b_vv1 = vv[row_next];
                    double b_ww1 = ww[row_next];
                    double b_t1  = time[row_next];
                    double b_duu = b_uu1 - b_uu;
                    double b_dvv = b_vv1 - b_vv;
                    double b_dww = b_ww1 - b_ww;
                    b_duvw = sqrt(b_duu * b_duu + b_dvv * b_dvv +
                        b_dww * b_dww);
                    b_dt = (b_t1 - b_t);
                }

                // If last time or if next point extends beyond average save out
                // average baseline and reset average, else accumulate current
                // average length.
                if (t == num_times - 1 || duvw[b] + b_duvw >= duvw_max
                                || dt[b] + b_dt >= dt_max) {
                    ms_out.addRow(1, false);

                    uInt bda_count = ave_count[b];

                    double bda_t  = ave_t[b] / bda_count;

                    // Populate coordinate column.
                    bda_uvw[0] = ave_uu[b] / bda_count;
                    bda_uvw[1] = ave_vv[b] / bda_count;
                    bda_uvw[2] = ave_ww[b] / bda_count;
                    col_uvw_out.put(bda_row, bda_uvw);
                    // Populate antenna columns.
                    col_ant1_out.put(bda_row, static_cast<int>(a1));
                    col_ant2_out.put(bda_row, static_cast<int>(a2));
                    // Populate flag column.
                    col_flag_out.put(bda_row, flag);
                    // Populate time columns.
                    col_time_out.put(bda_row, bda_t);
                    col_time_centroid_out.put(bda_row, bda_t);
                    col_interval_out.put(bda_row, bda_count * delta_t);
                    col_exposure_out.put(bda_row, bda_count * delta_t);
                    // Populate weight and sigma columns
                    float bda_sigma = 1.0 / sqrt(static_cast<float>(bda_count));
                    float bda_weight = static_cast<float>(bda_count);
                    col_sigma_out.put(bda_row, Vector<Float>(num_pols,
                        bda_sigma));
                    col_weight_out.put(bda_row, Vector<Float>(num_pols,
                        bda_weight));
                    // Populate data column(s).
                    for (int p = 0; p < num_pols; ++p) {
                        double bda_re = ave_data[b * num_pols + p].real() /
                            bda_count;
                        double bda_im = ave_data[b * num_pols + p].imag() /
                            bda_count;
                        bda_data_p[p] = Complex(bda_re, bda_im);
                    }
                    col_data_out.put(bda_row, bda_data);

                    // Reset baseline accumulation buffers.
                    duvw[b] = 0.0;
                    dt[b] = 0.0;
                    ave_count[b] = 0;
                    ave_uu[b] = 0.0;
                    ave_vv[b] = 0.0;
                    ave_ww[b] = 0.0;
                    ave_t[b] = 0.0;
                    for (int p = 0; p < num_pols; ++p) {
                        ave_data[b*num_pols + p].real() = 0.0;
                        ave_data[b*num_pols + p].imag() = 0.0;
                    }
                    // Update baseline average row counter.
                    bda_row += 1;
                }
                else {
                    // Accumulate distance and time on current baseline.
                    duvw[b] += b_duvw;
                    dt[b] += b_dt;
                }
            }
        }
    }
    gettimeofday(&timer_loop[1], 0);
    cout << "+ Loop timer elapsed   : " <<
        timer_elapsed(timer_loop[0], timer_loop[1]) << " s" << endl;

    gettimeofday(&timer_all[1], NULL);
    cout << "+ Total time elapsed   : " <<
        timer_elapsed(timer_all[0], timer_all[1]) << " s" << endl;

    cout << "+ Number of BDA rows   : " << ms_out.nrow() << endl;
    cout << "+ Data reduction       : " <<
        ms.nrow() / static_cast<double>(ms_out.nrow()) << ":1" << endl;

    cout << endl;
    cout << "BDA MS                 : " << filename_out << endl;

    return 0;
}
