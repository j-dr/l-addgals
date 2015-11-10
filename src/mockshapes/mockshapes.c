/*!\file mockshapes.c
 * \brief Add shapes and sizes to galaxy catalog based on lensfit/SuprimeCam models
 */
#include "mockshapes.h"

fitsfile *prepare_table(int *, long *, int *, long *);
int write_to_table(fitsfile *, double *, double *, int, long, long, long);
long read_magnitudes(fitsfile *, long, long, int, long, long, double *);
void add_size_shape();
int col_exists(fitsfile *fptr, char *colname);

int
main(int argc, char **argv)
{
    int a, narg, opt;
    static char prefsname[PATH_MAX];
    char **argkey, **argval;

    msg_init(MSG_INFO, NULL);
    if (argc < 2) {
        msg_message(MSG_INFO, 0, "%s-%s", BANNER, VERSION);
        msg_message(MSG_INFO, 0, "Written by %s", COPYRIGHT);
        msg_message(MSG_FATAL, 0, "%s", "SYNTAX:");
        msg_message(MSG_FATAL, 0, "%s", SYNTAX);
        exit(EXIT_FAILURE);
    }

    CALLOC(argkey, argc, "Could not allocate argkey: ");
    CALLOC(argval, argc, "Could not allocate argval: ");

    /*
     * default config file 
     */
    strcpy(prefsname, "mockshapes_default.conf");

    /*
     * This is all preferences and option reading taking from SWarp 
     */
    narg = 0;
    for(a = 1; a < argc; a++) {
        if (*(argv[a]) == '-') {
            opt = (int) argv[a][1];
            if (strlen(argv[a]) < 3 || opt == '-') {
                if (opt == '-')
                    opt = (int) tolower((int) argv[a][2]);
                switch (opt) {
                case 'c':
                    if (a < (argc - 1))
                        strcpy(prefsname, argv[++a]);
                    break;
                case 'd':
                    dumpprefs();
                    exit(EXIT_SUCCESS);
                    break;
                case 'h':
                default:
                    msg_message(MSG_FATAL, 0, "SYNTAX:\n%s", SYNTAX);
                }
            } else {
                argkey[narg] = &argv[a][1];
                argval[narg++] = argv[++a];
            }
        } else {
            snprintf(prefs.incat, PATH_MAX, "%s", argv[a]);
        }
    }

    readprefs(prefsname, argkey, argval, narg);

    msg_init(prefs.verbose, NULL);
    msg_message(MSG_INFO, 0, "%s-%s", BANNER, VERSION);
    msg_message(MSG_INFO, 0, "by %s", COPYRIGHT);
    msg_message(MSG_DEBUG, 0, "Verbosity level: %d", prefs.verbose);
#ifdef _OPENMP
    if (prefs.nthreads)
        omp_set_num_threads(prefs.nthreads);
    prefs.nthreads = omp_get_max_threads();
    msg_message(MSG_NOTICE, 0, "Running OpenMP version with %d threads.",
                prefs.nthreads);
#else
    msg_message(MSG_NOTICE, 0,
                "This build of mockshapes is single-threaded.");
    prefs.nthreads = 1;
#endif
    gsl_set_error_handler_off();
    add_size_shape();

    msg_finish();
    return EXIT_SUCCESS;
}



/*!\brief Open fits table, read table information
 * \param ncols total number of columns in table
 * \param nrows total number of rows in table
 * \param colnum column number of the magnitude column
 * \param vl length of magnitude vector
 * \return fits pointer to the input catalog, NULL in case of error
 */
fitsfile *
prepare_table(int *ncols, long *nrows, int *colnum, long *vl)
{
    int status = 0;
    int exist_colnum;
    int type;
    long width;
    char *ttype[] = { prefs.ename, prefs.sname };
    char *tform[] = { "2E", "E" };
    char fits_err_msg[FITS_ERR_LEN];
    fitsfile *fptr;

    fits_open_table(&fptr, prefs.incat, READWRITE, &status);
    if (status) {
        msg_message(MSG_FATAL, 0, "Could not open input table %s",
                    prefs.incat);
        fits_report_error(stderr, status);
        return NULL;
    }

    if ((exist_colnum = col_exists(fptr, prefs.ename)) != 0) {
	if (exist_colnum < 0)
	    return NULL;
	else
	    fits_delete_col(fptr, exist_colnum, &status);
    }
    if ((exist_colnum = col_exists(fptr, prefs.sname)) != 0) {
	if (exist_colnum < 0)
	    return NULL;
	else
	    fits_delete_col(fptr, exist_colnum, &status);
    }
    fits_get_num_rows(fptr, nrows, &status);
    fits_get_num_cols(fptr, ncols, &status);
    fits_get_colnum(fptr, 1, prefs.magname, colnum, &status);
    fits_get_coltype(fptr, *colnum, &type, vl, &width, &status);
    if (status) {
        fits_close_file(fptr, &status);
        msg_message(MSG_FATAL, 0, "Error accesscing table:");
        while (fits_read_errmsg(fits_err_msg))
            msg_message(MSG_FATAL, 0, fits_err_msg);
        return NULL;
    }
    msg_message(MSG_DEBUG, 0, "Table with %ld columns", *ncols);
    msg_message(MSG_DEBUG, 0, "Magnitude is in column %d", *colnum);
    msg_message(MSG_INFO, 0,
                "Table with %d rows and %d element magnitude vector",
                *nrows, *vl);

    fits_insert_cols(fptr, *ncols + 1, 2, ttype, tform, &status);
    if (status) {
	msg_message(MSG_FATAL, 0, "CFITSIO error code %d", status);
        msg_message(MSG_FATAL, 0, "Error adding columns to table:");
        fits_close_file(fptr, &status);
        while (fits_read_errmsg(fits_err_msg))
            msg_message(MSG_FATAL, 0, fits_err_msg);
        return NULL;
    }
    return fptr;
}


/*!\brief Write ellipticity and size information to fits table
 * \param fptr a pointer to an open fits table
 * \param e an array of length 2*blocksize holding pairs of e1, e2 components
 * \param s an array of length blocksize holding galaxy size in arcsec
 * \param ncols number of columns in the original table (new columns are appended)
 * \param iter count of the iterations
 * \param blocksize number of rows written/read at a time
 * \param nrows number of rows to be written in this call
 * \return 1 in case of success, 0 in case of error
 */
int
write_to_table(fitsfile * fptr, double *e, double *s, int ncols, long iter,
               long blocksize, long nrows)
{
    int status = 0;
    char fits_err_msg[FITS_ERR_LEN];

    fits_write_col(fptr, TDOUBLE, ncols + 1, iter * blocksize + 1, 1,
                   2 * nrows, e, &status);
    fits_write_col(fptr, TDOUBLE, ncols + 2, iter * blocksize + 1, 1, nrows,
                   s, &status);
    if (status) {
        fits_close_file(fptr, &status);
        msg_message(MSG_FATAL, 0, "Error writing to new columns:");
        while (fits_read_errmsg(fits_err_msg))
            msg_message(MSG_FATAL, 0, fits_err_msg);
        return 0;
    }
    return 1;
}


/*!\brief Read magnitude from fits table and return number of read objects 
 * \param fptr a pointer to an open fits table
 * \param vl length of the magnitude vector 
 * \param iter count of the iterations
 * \param colnum number of the column containing the magnitude vector
 * \param blocksize number of rows to be read in each iteration
 * \param nrows total number of rows in the table
 * \param mag pointer to an array of length blocksize into which the magnitudes will be written
 * \return number of objects read
 */
long
read_magnitudes(fitsfile * fptr, long vl, long iter, int colnum,
                long blocksize, long nrows, double *mag)
{
    int status = 0;
    long toread;

    if (iter * blocksize > nrows)
        return 0;
    if ((iter + 1) * blocksize > nrows)
        toread = nrows - iter * blocksize;
    else
        toread = blocksize;
    fits_read_col(fptr, TDOUBLE, colnum, iter * blocksize + 1, 1, vl * toread,
                  NULL, mag, NULL, &status);
    if (status) {
        msg_message(MSG_FATAL, 0, "Could not read magnitude values:");
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        return 0;
    }
    return toread;
}


/*!\brief main function to add sizes and shapes. Here is the iterator loop over the catalog.
 * \param None takes all input from the global prefs structure
 * \return Nothing. This function terminates the program in case of problems.
 */
void
add_size_shape()
{
    int status = 0;
    int ncols, colnum;
    int iam;
    long i, iter, nrows, nread, blocksize;
    long vl;                /* length of vector */
    double *mag, *e, *s;
    double mymag, e1, e2;
    eparam eparams;
    sparam sparams;
    gsl_rng **myrng;
    fitsfile *fptr;

    if (!(fptr = prepare_table(&ncols, &nrows, &colnum, &vl)))
        exit(EXIT_FAILURE);
    CALLOC(myrng, prefs.nthreads, "Random number generators");
    for(i = 0; i < prefs.nthreads; i++) {
        myrng[i] = gsl_rng_alloc(gsl_rng_ranlxd2);
        gsl_rng_set(myrng[i], prefs.seed + i);
    }
    CALLOC(e, 2 * nrows, "Ellipticity vector");
    CALLOC(s, nrows, "Size vector");

    /*
     * Determine the optimal number of rows to read/write at a time 
     */
    fits_get_rowsize(fptr, &blocksize, &status);
    msg_message(MSG_DEBUG, 0, "Will read chunks of %ld rows", blocksize);
    if (blocksize < prefs.nthreads) {
        msg_message(MSG_WARNING, 0,
                    "Line number (%d) is smaller than number of threads (%d)",
                    blocksize, prefs.nthreads);
        msg_message(MSG_WARNING, 0,
                    "Program will continue to run but inefficiently");
    }
    CALLOC(mag, vl * blocksize, "magnitude array");

    iter = 0;
    while ((nread = read_magnitudes(fptr, vl, iter, colnum, blocksize, nrows,
                                    mag))) {
        msg_message(MSG_INFO, 1, "Processing chunk %ld/%ld", iter + 1,
                    nrows / blocksize + 1);

#pragma omp parallel for default(none) shared(nrows, prefs, mag, vl, e, s, myrng, nread) private(iam, eparams, sparams, mymag, e1, e2)
        for(i = 0; i < nread; i++) {
#ifdef _OPENMP
            iam = omp_get_thread_num();
#else
            iam = 0;
#endif
            mymag = mag[prefs.nelem - 1 + i * vl];
            calceparams(mymag, &eparams);
            calcsparams(mymag, &sparams);
            rng_efunc(myrng[iam], &eparams, &e1, &e2);
            /*
             * These writes to shared arrays probably destroy cache
             * coherence
             */
            e[2 * i] = e1;
            e[2 * i + 1] = e2;
            s[i] = ran_gno(myrng[iam], &sparams);
            if (iam == 0 && !(i % 100))
                msg_message(MSG_DEBUG, 1, "Computed %d/%d", i,
                            nread / prefs.nthreads);
        }
        if (!write_to_table(fptr, e, s, ncols, iter, blocksize, nread))
            exit(EXIT_FAILURE);
        iter++;
        msg_message(MSG_DEBUG, 2, "Done");
    }
    msg_message(MSG_INFO, 2, "Done");
    for(i = 0; i < prefs.nthreads; i++)
        gsl_rng_free(myrng[i]);

    fits_close_file(fptr, &status);
}
