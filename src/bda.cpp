/* Copyright (c) 2015, The University of Oxford
   BDS 3-Clause license (see LICENSE) */
#include <casacore/tables/Tables.h>
#include <casacore/ms/MeasurementSets.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <string>
#include <complex>
#include <cmath>
#include <vector>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using casa::Table;
using casa::TableColumn;
using casa::TableExprNode;
using casa::Array;
using casa::Double;
using casa::Vector;
using casa::uInt;
using casa::ArrayColumn;
using casa::MeasurementSet;
using casa::IPosition;
using casa::ROMSColumns;

using namespace casa;

/**
 * arcsinc(x) function taken from Obit. Uses Newton-Raphson method.
 */
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
 * - Handle polarisation
 * - Copy MS sub-tables (antenna etc)
 * - Evaluate stats. on averages in the same way as mstransform()
 */

int main(int argc, char** argv) {
    if (argc != 5) {
        cout << "Usage:" << endl;
        cout << "  ./bda <ms> <max_fact> <fov> <idt_max>" << endl << endl;
        cout << "Where:" << endl;
        cout << "  ms       Path of measurement set to be averaged." << endl;
        cout << "  max_fact Maximum amplitude reduction factor. (eg. 1.01)" << endl;
        cout << "  fov      Field-of-view radius in degrees associated " << endl;
        cout << "           with the amplitude reduction factor" << endl;
        cout << "  idt_max  Maximum averaging time in terms of non-averaged" << endl;
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

    struct timeval timer_all[2];
    gettimeofday(&timer_all[0], NULL);

    // const MeasurementSet ms(filename);
    const Table ms(filename);
    uInt num_rows = ms.nrow();

    // Make a copy of the table for BDA.
    bool noRows = true;
    bool valueCopy = false;
    ms.deepCopy(filename_out, Table::New, valueCopy, Table::LittleEndian,
                noRows);
    Table ms_out(filename_out, Table::Update);
    ArrayColumn<Double> col_uvw_out(ms_out, "UVW");
    ScalarColumn<Double> col_time_out(ms_out, "TIME");
    ScalarColumn<Double> col_time_centroid_out(ms_out, "TIME_CENTROID");
    ScalarColumn<Int> col_ant1_out(ms_out, "ANTENNA1");
    ScalarColumn<Int> col_ant2_out(ms_out, "ANTENNA2");
    ScalarColumn<Double> col_interval_out(ms_out, "INTERVAL");
    ScalarColumn<Double> col_exposure_out(ms_out, "EXPOSURE");
    ArrayColumn<Complex> col_data_out(ms_out, "DATA");
    ArrayColumn<Float> col_weight_out(ms_out, "WEIGHT");
    ArrayColumn<Float> col_sigma_out(ms_out, "SIGMA");


    // Load antenna table to get number of antennas.
    const Table ant(filename+"/ANTENNA");
    uInt num_antennas = ant.nrow();
    ant.deepCopy(filename_out+"/ANTENNA", Table::New, true,
                 Table::LittleEndian, false);

    // Work out number of baselines and number of times.
    uInt num_baselines = num_antennas * (num_antennas - 1) / 2;
    uInt num_times = num_rows / num_baselines;

    // Load observation table to get time range (needed for delta_t).
    const Table obs(filename+"/OBSERVATION");
    obs.deepCopy(filename_out+"/OBSERVATION", Table::New, true,
                 Table::LittleEndian, false);
    const ArrayColumn<Double> time_range(obs, "TIME_RANGE");
    const Array<Double> time_range_cell = time_range.get(0);
    double t0 = time_range_cell.data()[0];
    double t1 = time_range_cell.data()[1];
    double obs_length_s = t1 - t0;
    double delta_t = obs_length_s / num_times;

    // Load frequency from MS. (require 1 freq)
    const Table spectral_window(filename+"/SPECTRAL_WINDOW");
    spectral_window.deepCopy(filename_out+"/SPECTRAL_WINDOW", Table::New, true,
                             Table::LittleEndian, false);
    const ArrayColumn<Double> chan_freq(spectral_window, "CHAN_FREQ");
    const Array<Double> freqs0 = chan_freq.get(0);
    double freq = freqs0.data()[0];
    double wavelength = 299792458.0 / freq;

    // Copy other sub-tables. (ANTENNA, OBSERVATION, SPECTRAL_WINDOW already
    // done)
    vector<string> sub_tables;
    sub_tables.push_back("/DATA_DESCRIPTION");
    sub_tables.push_back("/FEED");
    sub_tables.push_back("/FIELD");
    sub_tables.push_back("/HISTORY");
    sub_tables.push_back("/POLARIZATION");
    for (vector<string>::iterator it = sub_tables.begin();
                    it != sub_tables.end(); ++it) {
        const Table sub_table(filename+*it);
        sub_table.deepCopy(filename_out+*it, Table::New, true,
                           Table::LittleEndian, false);
    }

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
    printf("\n");

    const ArrayColumn<Double> col_uvw(ms, "UVW");
    const IPosition uvw_shape = col_uvw.shapeColumn();
    const ScalarColumn<Double> col_time(ms, "TIME");
    const ArrayColumn<Complex> col_data(ms, "DATA");
    const IPosition data_shape = col_data.shapeColumn();

    std::vector<double> uu_(num_rows);
    std::vector<double> vv_(num_rows);
    std::vector<double> ww_(num_rows);
    std::vector<double> time_(num_rows);
    std::vector<Complex> data_(num_rows);
    double* uu = &uu_[0];
    double* vv = &vv_[0];
    double* ww = &ww_[0];
    double* time = &time_[0];
    Complex* data = &data_[0];

    struct timeval timer_load[2];
    gettimeofday(&timer_load[0], NULL);
    Array<Double> temp_uvw(uvw_shape);
    // TODO(BM) deal with polarised data files.
    for (uInt i = 0; i < num_rows; ++i) {
        temp_uvw = col_uvw.get(i);
        uu[i] = temp_uvw.data()[0];
        vv[i] = temp_uvw.data()[1];
        ww[i] = temp_uvw.data()[2];
        time[i] = col_time.get(i);
        data[i] = col_data.get(i).data()[0];
    }
    gettimeofday(&timer_load[1], NULL);
    double elapsed = (timer_load[1].tv_sec - timer_load[0].tv_sec) +
                    (timer_load[1].tv_usec - timer_load[0].tv_usec) / 1.0e6;
    cout << "+ Loading data took    : " << elapsed << " s" << endl;

    // Buffers of deltas along the baseline in the current average.
    std::vector<double> duvw_(num_baselines, 0.0);
    std::vector<double> dt_(num_baselines, 0);
    double* duvw = &duvw_[0];
    double* dt = &dt_[0];
    // Buffers of the current average along the baseline.
    std::vector<double> ave_uu_(num_baselines, 0.0);
    std::vector<double> ave_vv_(num_baselines, 0.0);
    std::vector<double> ave_ww_(num_baselines, 0.0);
    std::vector<double> ave_t_(num_baselines, 0.0);
    std::vector<Complex> ave_data_(num_baselines, 0.0);
    std::vector<uInt> ave_count_(num_baselines, 0);
    double* ave_uu = &ave_uu_[0];
    double* ave_vv = &ave_vv_[0];
    double* ave_ww = &ave_ww_[0];
    double* ave_t = &ave_t_[0];
    uInt* ave_count = &ave_count_[0];
    Array<Double> bda_uvw(uvw_shape);
    Array<Complex> bda_data(data_shape);
    Complex* ave_data = &ave_data_[0];
    uInt num_pols = 1;  // TODO(BM) get this from the MS and use elsewhere too.
    Vector<Float> bda_weight(num_pols);

    struct timeval timer_loop[2];
    gettimeofday(&timer_loop[0], 0);
    for (uInt t = 0; t < num_times; ++t) {
        // for (uInt b = 0; b < num_baselines; ++b) {
        for (uInt b = 0, a1 = 0; a1 < num_antennas; ++a1) {
            for (uInt a2 = a1+1; a2 < num_antennas; ++a2, ++b) {

                uInt row = t * num_baselines + b;

                // Get coordinates and time for the current baseline.
                double b_uu = uu[row];
                double b_vv = vv[row];
                double b_ww = ww[row];
                double b_t  = time[row];
                Complex b_data = data[row];

                // Accumulate into the average.
                ave_uu[b] += b_uu;
                ave_vv[b] += b_vv;
                ave_ww[b] += b_ww;
                ave_t[b] += b_t;
                ave_data[b] += b_data;
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
                    b_duvw = sqrt(b_duu * b_duu + b_dvv * b_dvv + b_dww * b_dww);
                    b_dt = (b_t1 - b_t);
                }

//                if (b == 0) {
//                    printf("  b:%i t:%i b_duvw:%f duvw:%f b_dt:%f dt:%f, b_uu:%f\n",
//                           b, t, b_duvw, duvw[b], b_dt, dt[b], b_uu);
//                }

                // If last time or if next point extends beyond average save out
                // average baseline and reset average, else accumulate current
                // average length.
                if (t == num_times - 1 || duvw[b]+b_duvw >= duvw_max
                                || dt[b] + b_dt >= dt_max) {

                    uInt bda_count = ave_count[b];
                    double bda_uu = ave_uu[b] / bda_count;
                    double bda_vv = ave_vv[b] / bda_count;
                    double bda_ww = ave_ww[b] / bda_count;
                    double bda_t  = ave_t[b] / bda_count;
                    double bda_re = ave_data[b].real() / bda_count;
                    double bda_im = ave_data[b].imag() / bda_count;
                    bda_data[0] = Complex(bda_re, bda_im);
                    ms_out.addRow(1, true);
                    uInt outrow = ms_out.nrow() - 1;
                    bda_uvw[0] = bda_uu;
                    bda_uvw[1] = bda_vv;
                    bda_uvw[2] = bda_ww;
                    col_uvw_out.put(outrow, bda_uvw);
                    col_ant1_out.put(outrow, (int)a1);
                    col_ant2_out.put(outrow, (int)a2);
                    col_time_out.put(outrow, bda_t);
                    col_time_centroid_out.put(outrow, bda_t);
                    col_data_out.put(outrow, bda_data);
                    col_interval_out.put(outrow, bda_count * delta_t);
                    col_exposure_out.put(outrow, bda_count * delta_t);
                    bda_weight[0] = (float)ave_count[b];
                    col_weight_out.put(outrow, bda_weight);
                    double sigma = 1.0 / sqrt(static_cast<float>(bda_count));
                    col_sigma_out.put(outrow, Vector<Float>(num_pols, sigma));

//                    if (b == 0) {
//                        printf("  count:%i uu:(%f) %f\n", ave_count[b], bda_uu, ave_uu[b]);
//                        cout << endl;
//                    }
                    duvw[b] = 0.0;
                    dt[b] = 0.0;
                    ave_count[b] = 0;
                    ave_uu[b] = 0.0;
                    ave_vv[b] = 0.0;
                    ave_ww[b] = 0.0;
                    ave_data[b] = Complex(0.0, 0.0);
                }
                else {
                    duvw[b] += b_duvw;
                    dt[b] += b_dt;
                }
            }
        }
    }
    gettimeofday(&timer_loop[1], 0);
    elapsed = (timer_loop[1].tv_sec - timer_loop[0].tv_sec) +
                    (timer_loop[1].tv_usec - timer_loop[0].tv_usec) / 1.0e6;
    cout << "+ Loop timer elapsed   : " << elapsed << " s" << endl;
    cout << "  - Number of BDA rows : " << ms_out.nrow() << endl;
    cout << "  - Data reduction     : " << ms.nrow()/(double)ms_out.nrow() << ":1" << endl;

    gettimeofday(&timer_all[1], NULL);
    elapsed = (timer_all[1].tv_sec - timer_all[0].tv_sec) +
                    (timer_all[1].tv_usec - timer_all[0].tv_usec) / 1.0e6;
    cout << "+ Total time elapsed   : " << elapsed << " s" << endl;

    return 0;
}
