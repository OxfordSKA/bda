/* Copyright (c) 2015, The University of Oxford
   BDS 3-Clause license (see LICENSE) */
#include <casacore/tables/Tables.h>
#include <casacore/ms/MeasurementSets.h>
#include <cstdio>
#include <iostream>
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

int main(int argc, char** argv) {
    if (argc != 2) {
        cout << "Usage: bda <ms>" << endl;
        return 0;
    }
    string filename = string(argv[1]);
    string filename_out = filename+".bda";

#if 0

    Table t(filename);
    TableColumn col_uvw = TableColumn(t, "UVW");
    TableExprNode col_uvw_ = t.col("UVW");

    cout << "No. rows  : " << col_uvw.nrow() << endl;
    ArrayColumn<Double> a_uvw(t, "UVW");

    const Array<Double> x = a_uvw.get(0);
    double uu = *(x[0].data());
    double vv = *(x[1].data());
    double ww = *(x[2].data());
    cout << uu << " " << vv << " " << ww << endl;

    // TableColumn col_time = TableColumn(t, "TIME");

#else

    struct timeval timer_all[2];
    gettimeofday(&timer_all[0], NULL);

    // MeasurementSet ms(filename);
    const Table ms(filename);
    uInt num_rows = ms.nrow();

    // Make a copy of the table for BDA.
    bool noRows = true;
    bool valueCopy = false;
    ms.deepCopy(filename_out, Table::New, valueCopy, Table::LittleEndian,
                noRows);
    // TODO(BM) copy sub-tables (antenna table etc.).
    Table ms_out(filename_out, Table::Update);
    ArrayColumn<Double> col_uvw_out(ms_out, "UVW");
    ScalarColumn<Double> col_time_out(ms_out, "TIME");
    ScalarColumn<Int> col_ant1_out(ms_out, "ANTENNA1");
    ScalarColumn<Int> col_ant2_out(ms_out, "ANTENNA2");
    ArrayColumn<Complex> col_data_out(ms_out, "DATA");

    // Load antenna table to get number of antennas.
    const Table ant(filename+"/ANTENNA");
    uInt num_antennas = ant.nrow();

    // Work out number of baselines and number of times.
    uInt num_baselines = num_antennas * (num_antennas - 1) / 2;
    uInt num_times = num_rows / num_baselines;

    // Load observation table to get time range (needed for delta_t).
    Table obs(filename+"/OBSERVATION");
    ArrayColumn<Double> time_range(obs, "TIME_RANGE");
    double t0 = *static_cast<double*>(time_range.get(0)[0].data());
    double t1 = *static_cast<double*>(time_range.get(0)[1].data());
    double obs_length_s = t1 - t0;
    double delta_t = obs_length_s / num_times;

    printf("%s\n", string(80, '-').c_str());
    printf("MS            : %s\n", filename.c_str());
    printf("No. rows      : %i\n", num_rows);
    printf("No. antennas  : %i\n", num_antennas);
    printf("No. baselines : %i\n", num_baselines);
    printf("No. times     : %i\n", num_times);
    printf("delta_t       : %.16f\n", delta_t);
    printf("%s\n", string(80, '-').c_str());

    const ArrayColumn<Double> col_uvw(ms, "UVW");
    const ScalarColumn<Double> col_time(ms, "TIME");
    const ArrayColumn<Complex> col_data(ms, "DATA");

    std::vector<double> uu_(num_rows);
    std::vector<double> vv_(num_rows);
    std::vector<double> ww_(num_rows);
    std::vector<double> time_(num_rows);
    std::vector<Complex> data_(num_rows);
    double* uu = &uu_[0];
    double* vv = &uu_[0];
    double* ww = &uu_[0];
    double* time = &time_[0];
    Complex* data = &data_[0];

    struct timeval timer_load[2];
    gettimeofday(&timer_load[0], NULL);
    for (uInt i = 0; i < num_rows; ++i) {
        uu[i] = *(col_uvw.get(i)[0].data());
        vv[i] = *(col_uvw.get(i)[1].data());
        ww[i] = *(col_uvw.get(i)[2].data());
        time[i] = col_time.get(i);
        data[i] = *(col_data.get(i).data());
    }
    gettimeofday(&timer_load[1], NULL);
    double elapsed = (timer_load[1].tv_sec - timer_load[0].tv_sec) +
                    (timer_load[1].tv_usec - timer_load[0].tv_usec) / 1.0e6;
    cout << "+ Loading data took  : " << elapsed << " s" << endl;

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

    Complex* ave_data = &ave_data_[0];

    // TODO(BM) work these out properly..
    double dt_max = 4 * delta_t;
    double duvw_max = 2.2;

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

                if (b == 0) {
                    printf("  b:%i t:%i b_duvw:%f duvw:%f b_dt:%f dt:%f, b_uu:%f\n",
                           b, t, b_duvw, duvw[b], b_dt, dt[b], b_uu);
                }

                // If last time or if next point extends beyond average save out
                // average baseline and reset average, else accumulate current
                // average length.
                if (t == num_times - 1 || duvw[b]+b_duvw >= duvw_max
                                || dt[b] + b_dt >= dt_max) {

                    double bda_uu = ave_uu[b] / ave_count[b];
                    double bda_vv = ave_vv[b] / ave_count[b];
                    double bda_ww = ave_ww[b] / ave_count[b];
                    double bda_t  = ave_t[b] / ave_count[b];
                    double bda_re = ave_data[b].real() / ave_count[b];
                    double bda_im = ave_data[b].imag() / ave_count[b];
                    Complex bda_data(bda_re, bda_im);
                    ms_out.addRow(1, true);
                    uInt outrow = ms_out.nrow() - 1;
                    col_time_out.put(outrow, bda_t);
                    col_ant1_out.put(outrow, a1);
                    col_ant2_out.put(outrow, a2);
                    if (b == 0) {
                        printf("  count:%i uu:(%f) %f\n", ave_count[b], bda_uu, ave_uu[b]);
                        cout << endl;
                    }

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
    cout << "+ Loop timer elapsed : " << elapsed << " s" << endl;
    cout << "+ Number of BDA rows : " << ms_out.nrow() << endl;
    cout << "  - Data reduction   : " << ms.nrow()/(double)ms_out.nrow() << ":1" << endl;

    gettimeofday(&timer_all[1], NULL);
    elapsed = (timer_all[1].tv_sec - timer_all[0].tv_sec) +
                    (timer_all[1].tv_usec - timer_all[0].tv_usec) / 1.0e6;
    cout << "+ Total time elapsed : " << elapsed << " s" << endl;

#endif

    return 0;
}
