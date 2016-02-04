/* Copyright (c) 2015-2016, The University of Oxford
   BDS 3-Clause license (see LICENSE) */
#include <sys/time.h>
#include <casacore/tables/Tables.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <cassert>

using namespace casa;

double timer_elapsed(const struct timeval& t1, const struct timeval& t2)
{
    return (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1.0e6;
}

void copy_ms(const string& ms_in, const string& ms_out)
{
    // Copy the main table headers.
    const Table ms(ms_in);
    ms.deepCopy(ms_out, Table::New, false, Table::LittleEndian, true);

    // Copy all the sub-tables with their data, except SORTED_TABLE.
    uInt num_fields = ms.keywordSet().description().nfields();
    for (uInt i = 0; i < num_fields; ++i)
    {
        if (ms.keywordSet().description().isTable(i))
        {
            string table_name = ms.keywordSet().description().name(i);
            if (table_name == "SORTED_TABLE") continue;
            const Table sub_table(ms_in + "/" + table_name);
            sub_table.deepCopy(ms_out + "/" + table_name, Table::New, true,
                               Table::LittleEndian, false);
        }
    }
}


/* arcsinc(x) function taken from Obit. Uses Newton-Raphson method.  */
double inv_sinc(double value)
{
    double x1 = 0.001;
    for (int i = 0; i < 1000; ++i)
    {
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
 * - Exploit extra performance of blocks and addrow(init=false)
 */
int main(int argc, char** argv)
{
    if (argc != 6)
    {
        printf("Usage:\n");
        printf("  ./bda <ms_in> <ms_out> <max_fact> <fov> <dt_max>\n\n");
        printf("Where:\n");
        printf("  ms_in    Path of measurement set to be averaged.\n");
        printf("  ms_out   Path of output measurement set.\n");
        printf("  max_fact Maximum amplitude reduction factor. (eg. 1.01)\n");
        printf("  fov      Field-of-view radius in degrees associated\n");
        printf("           with the amplitude reduction factor.\n");
        printf("  dt_max   Maximum averaging time in seconds.\n");
        printf("Example:\n");
        printf("  ./bda vis.ms vis_bda.ms 1.01 0.9 10.0\n");
        return 0;
    }
    string filename = string(argv[1]);
    string filename_out = string(argv[2]);
    double max_fact = atof(argv[3]);
    double fov = atof(argv[4]);
    double dt_max = atof(argv[5]);

    // Create a total run time timer.
    struct timeval tmr[2];
    gettimeofday(&tmr[0], NULL);

    // Copy input measurement set.
    copy_ms(filename, filename_out);
    const Table ms(filename);
    uInt num_rows = ms.nrow();
    if (num_rows == 0)
    {
        fprintf(stderr, "ERROR: No data in input Measurement Set!\n");
        return EXIT_FAILURE;
    }

    // Get number of antennas, number of baselines and number of times.
    uInt num_antennas = 0;
    {
        const Table ant(filename + "/ANTENNA");
        num_antennas = ant.nrow();
    }
    uInt num_baselines = num_antennas * (num_antennas - 1) / 2;
    uInt num_times = num_rows / num_baselines;

    // Get delta_t.
    double delta_t = 0.0;
    {
        const Table t(filename + "/OBSERVATION");
        const ArrayColumn<Double> time_range(t, "TIME_RANGE");
        const Array<Double> cell = time_range.get(0);
        double obs_length_s = cell.data()[1] - cell.data()[0];
        delta_t = obs_length_s / num_times;
    }

    // Get frequency.
    double freq = 0.0;
    {
        const Table t(filename + "/SPECTRAL_WINDOW");
        const ArrayColumn<Double> col_chan_freq(t, "CHAN_FREQ");
        const Vector<Double> freqs0 = col_chan_freq.get(0);
        int num_freqs = freqs0.size();
        if (num_freqs != 1) {
            fprintf(stderr, "ERROR: BDA only works with MS containing "
                "one frequency! (Found %i frequencies)\n", num_freqs);
            return EXIT_FAILURE;
        }
        freq = freqs0.data()[0];
    }
    double wavelength = 299792458.0 / freq;

    // Get number of polarisations.
    int num_pols = 0;
    {
        const Table t(filename + "/POLARIZATION");
        const ArrayColumn<Int> corr_type(t, "CORR_TYPE");
        num_pols = corr_type.get(0).size();
    }

    // Create references to standard columns.
    const ArrayColumn<Double> col_uvw(ms, "UVW");
    const ArrayColumn<Complex> col_data(ms, "DATA");
    const ScalarColumn<Double> col_time(ms, "TIME");

    // Check if extra columns are present.
    ArrayColumn<Complex> col_model_data, col_corrected_data;
    const bool have_model_data = ms.tableDesc().isColumn("MODEL_DATA");
    const bool have_corrected_data = ms.tableDesc().isColumn("CORRECTED_DATA");
    if (have_model_data)
        col_model_data.attach(ms, "MODEL_DATA");
    if (have_corrected_data)
        col_corrected_data.attach(ms, "CORRECTED_DATA");

    // ------------------------------------------------------------------------

    // Evaluate averaging parameters.
    double duvw_max = inv_sinc(1.0 / max_fact) / (fov * (M_PI / 180.0));
    duvw_max *= wavelength;

    // Print summary.
    printf(" | + Input MS              : %s\n", filename.c_str());
    printf(" |   + Number of rows      : %i\n", num_rows);
    printf(" |   + Number of antennas  : %i\n", num_antennas);
    printf(" |   + Number of baselines : %i\n", num_baselines);
    printf(" |   + Number of times     : %i\n", num_times);
    printf(" |   + Number of pols      : %i\n", num_pols);
    printf(" |   + Frequency           : %.6f MHz\n", freq/1.e6);
    printf(" |   + Time interval       : %.4f s\n", delta_t);
    printf(" | + Output MS             : %s\n", filename_out.c_str());
    printf(" |   + max_fact            : %f\n", max_fact);
    printf(" |   + Field of view       : %f deg\n", fov);
    printf(" |   + Maximum time        : %.4f s\n", dt_max);
    printf(" |   + Maximum UVW distance: %f m\n", duvw_max);
    printf(" | %s\n", string(77, '-').c_str());

    // Open the output MS and get references to output columns.
    Table ms_out(filename_out, Table::Update);
    ArrayColumn<Double> col_uvw_out(ms_out, "UVW");
    ScalarColumn<Double> col_time_out(ms_out, "TIME");
    ScalarColumn<Double> col_time_centroid_out(ms_out, "TIME_CENTROID");
    ScalarColumn<Int> col_ant1_out(ms_out, "ANTENNA1");
    ScalarColumn<Int> col_ant2_out(ms_out, "ANTENNA2");
    ArrayColumn<Bool> col_flag_out(ms_out, "FLAG");
    ScalarColumn<Double> col_interval_out(ms_out, "INTERVAL");
    ScalarColumn<Double> col_exposure_out(ms_out, "EXPOSURE");
    ArrayColumn<Float> col_weight_out(ms_out, "WEIGHT");
    ArrayColumn<Float> col_sigma_out(ms_out, "SIGMA");
    ArrayColumn<Complex> col_data_out(ms_out, "DATA");
    ArrayColumn<Complex> col_model_data_out, col_corrected_data_out;
    if (have_model_data)
        col_model_data_out.attach(ms_out, "MODEL_DATA");
    if (have_corrected_data)
        col_corrected_data_out.attach(ms_out, "CORRECTED_DATA");

    // Buffers of current and next UVW coordinates and times.
    IPosition len(1, num_baselines);
    IPosition uvw_size(2, 3, num_baselines);
    IPosition data_size(3, num_pols, 1, num_baselines);
    Array<Double> uvw_current_(uvw_size);
    Array<Double> uvw_next_(uvw_size);
    Vector<Double> time_current_(len);
    Vector<Double> time_next_(len);
    double* uvw_current = uvw_current_.data();
    double* uvw_next = uvw_next_.data();
    double* time_current = time_current_.data();
    double* time_next = time_next_.data();

    // Buffers of current visibility data.
    Array<Complex> data_(data_size);
    Array<Complex> model_data_(data_size);
    Array<Complex> corrected_data_(data_size);
    const Complex* data = data_.data();
    const Complex* model_data = model_data_.data();
    const Complex* corrected_data = corrected_data_.data();

    // Buffers of deltas along the baseline in the current average.
    double* duvw         = (double*) calloc(num_baselines, sizeof(double));
    double* dt           = (double*) calloc(num_baselines, sizeof(double));

    // Buffers of the current average along the baseline.
    double* ave_uu       = (double*) calloc(num_baselines, sizeof(double));
    double* ave_vv       = (double*) calloc(num_baselines, sizeof(double));
    double* ave_ww       = (double*) calloc(num_baselines, sizeof(double));
    double* ave_t        = (double*) calloc(num_baselines, sizeof(double));
    uInt* ave_count      = (uInt*)   calloc(num_baselines, sizeof(uInt));
    Complex* ave_data = 
            (Complex*) calloc(num_baselines * num_pols, sizeof(Complex));
    Complex* ave_model_data = 
            (Complex*) calloc(num_baselines * num_pols, sizeof(Complex));
    Complex* ave_corrected_data = 
            (Complex*) calloc(num_baselines * num_pols, sizeof(Complex));

    // CASA containers to hold one average only, for writing to output columns.
    Vector<Double> bda_uvw(3);
    Matrix<Complex> bda_data(num_pols, 1);
    Matrix<Complex> bda_model_data(num_pols, 1);
    Matrix<Complex> bda_corrected_data(num_pols, 1);
    Matrix<Bool> flag(num_pols, 1, false);
    Complex* bda_data_p = bda_data[0].data();
    Complex* bda_model_data_p = bda_model_data[0].data();
    Complex* bda_corrected_data_p = bda_corrected_data[0].data();

    // Read the first block of baseline coordinates.
    Slicer row_range(IPosition(1, 0), len);
    col_uvw.getColumnRange(row_range, uvw_current_);
    col_time.getColumnRange(row_range, time_current_);

    // Do the averaging.
    int bda_row = 0;
    for (uInt t = 0; t < num_times; ++t) 
    {
        // Read the visibility data for the current time.
        Slicer row_range(IPosition(1, t * num_baselines), len);
        col_data.getColumnRange(row_range, data_);
        if (have_model_data)
            col_model_data.getColumnRange(row_range, model_data_);
        if (have_corrected_data)
            col_corrected_data.getColumnRange(row_range, corrected_data_);

        // Make sure the column dimension matches the number of polarisations.
        assert(col_data.shapeColumn()[0] == num_pols);
        assert(col_data.shapeColumn()[1] == 1);

        // Read the next block of baseline coordinates if it exists.
        if (t < num_times - 1)
        {
            Slicer row_range(IPosition(1, (t + 1) * num_baselines), len);
            col_uvw.getColumnRange(row_range, uvw_next_);
            col_time.getColumnRange(row_range, time_next_);
        }

#if 0
        for (uInt i = 0; i < num_baselines; ++i)
        {
            printf("t = %2d, t_c = %.16f, t_n = %.16f, u_c = %10.3f, u_n = %10.3f, v_c = %10.3f, v_n = %10.3f, w_c = %10.3f, w_n = %10.3f\n",
                t, time_current[i], time_next[i], uvw_current[3*i + 0], uvw_next[3*i + 0], uvw_current[3*i + 1], uvw_next[3*i + 1], uvw_current[3*i + 2], uvw_next[3*i + 2]);
        }
#endif

        // Loop over baselines.
        for (uInt b = 0, a1 = 0; a1 < num_antennas; ++a1)
        {
            for (uInt a2 = a1 + 1; a2 < num_antennas; ++a2, ++b)
            {
                // Accumulate into averages.
                double b_uu = uvw_current[b*3 + 0];
                double b_vv = uvw_current[b*3 + 1];
                double b_ww = uvw_current[b*3 + 2];
                ave_count[b] += 1;
                ave_t[b]     += time_current[b];
                ave_uu[b]    += b_uu;
                ave_vv[b]    += b_vv;
                ave_ww[b]    += b_ww;
                for (int p = 0; p < num_pols; ++p)
                {
                    uInt i = b * num_pols + p;
                    ave_data[i] += data[i];
                }
                if (have_model_data)
                {
                    for (int p = 0; p < num_pols; ++p)
                    {
                        uInt i = b * num_pols + p;
                        ave_model_data[i] += model_data[i];
                    }
                }
                if (have_corrected_data)
                {
                    for (int p = 0; p < num_pols; ++p)
                    {
                        uInt i = b * num_pols + p;
                        ave_corrected_data[i] += corrected_data[i];
                    }
                }

                // Look ahead to the next point on the baseline and see if
                // this is also in the average.
                double b_duvw = 0.0, b_dt = 0.0;
                if (t < num_times - 1) 
                {
                    double du = uvw_next[b*3 + 0] - b_uu;
                    double dv = uvw_next[b*3 + 1] - b_vv;
                    double dw = uvw_next[b*3 + 2] - b_ww;
                    b_duvw = sqrt(du*du + dv*dv + dw*dw);
                    b_dt = time_next[b] - time_current[b];
                }

                // If last time or if next point extends beyond average,
                // save out averaged data and reset, 
                // else accumulate current average lengths.
                if (t == num_times - 1 || duvw[b] + b_duvw > duvw_max
                                || dt[b] + b_dt > dt_max) 
                {
                    ms_out.addRow(1, false);

                    // Populate meta-data columns.
                    uInt bda_count = ave_count[b];
                    double s = 1.0 / bda_count;
                    double bda_t  = ave_t[b] * s;
                    bda_uvw[0] = ave_uu[b] * s;
                    bda_uvw[1] = ave_vv[b] * s;
                    bda_uvw[2] = ave_ww[b] * s;
                    col_uvw_out.put(bda_row, bda_uvw);
                    col_ant1_out.put(bda_row, (int)a1);
                    col_ant2_out.put(bda_row, (int)a2);
                    col_flag_out.put(bda_row, flag);
                    col_time_out.put(bda_row, bda_t);
                    col_time_centroid_out.put(bda_row, bda_t);
                    col_interval_out.put(bda_row, bda_count * delta_t);
                    col_exposure_out.put(bda_row, bda_count * delta_t);
                    col_sigma_out.put(bda_row, 
                            Vector<Float>(num_pols, sqrt(s)));
                    col_weight_out.put(bda_row, 
                            Vector<Float>(num_pols, (float)bda_count));

                    // Populate data column(s).
                    for (int p = 0; p < num_pols; ++p) 
                    {
                        uInt i = b * num_pols + p;
                        bda_data_p[p] = ave_data[i] * s;
                        bda_model_data_p[p] = ave_model_data[i] * s;
                        bda_corrected_data_p[p] = ave_corrected_data[i] * s;
                    }
                    col_data_out.put(bda_row, bda_data);
                    if (have_model_data)
                        col_model_data_out.put(bda_row, bda_model_data);
                    if (have_corrected_data)
                        col_corrected_data_out.put(bda_row, bda_corrected_data);

                    // Reset baseline accumulation buffers.
                    duvw[b] = 0.0;
                    dt[b] = 0.0;
                    ave_count[b] = 0;
                    ave_t[b] = 0.0;
                    ave_uu[b] = 0.0;
                    ave_vv[b] = 0.0;
                    ave_ww[b] = 0.0;
                    for (int p = 0; p < num_pols; ++p) 
                    {
                        int i = b * num_pols + p;
                        ave_data[i] = Complex(0.0);
                        ave_model_data[i] = Complex(0.0);
                        ave_corrected_data[i] = Complex(0.0);
                    }

                    // Update baseline average row counter for next output.
                    ++bda_row;
                }
                else 
                {
                    // Accumulate distance and time on current baseline.
                    duvw[b] += b_duvw;
                    dt[b] += b_dt;
                }
            }
        }

        // Copy "next" to "current" for the next time block.
        memcpy(uvw_current, uvw_next, 3 * num_baselines * sizeof(double));
        memcpy(time_current, time_next, num_baselines * sizeof(double));
    }

    // Clean up.
    free(duvw);
    free(dt);
    free(ave_t);
    free(ave_uu);
    free(ave_vv);
    free(ave_ww);
    free(ave_count);
    free(ave_data);
    free(ave_model_data);
    free(ave_corrected_data);

    // Print summary.
    gettimeofday(&tmr[1], NULL);
    printf(" | + Elapsed time          : %.3f s\n", 
            timer_elapsed(tmr[0], tmr[1]));
    printf(" | + Number of BDA rows    : %i\n", ms_out.nrow());
    printf(" | + Data reduction        : %.2f:1\n",
            ms.nrow() / (double)ms_out.nrow());

    return EXIT_SUCCESS;
}
