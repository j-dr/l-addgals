#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <fitsio.h>

#include "define.h"
#include "message.h"
#include "prefs.h"


int col_exists(fitsfile *fptr, char *colname);


/*!\brief Add lensing columns to DES mock catalog
 * \param fname Filename of the DES mock catalog
 * \param nrows Pointer to number of rows in catalog
 * \return fits file pointer to mock catalog, NULL in case of error
 */
fitsfile *
prepare_table(char *fname, long *nrows)
{
    int i;
    int ncols, type, colnum;
    long width, repeat;
    int status = 0;
    /* Make sure that lensmagname is first. Other columns will be
       deleted for non-shape catalogs and deleting from the end is
       more efficient. */
    char *ttype[] = {prefs.lensmagname, "ORA", "ODEC", "GAMMA", "KAPPA", "MU", 
		     prefs.shearname, prefs.lenssname,};
    char *tform[] = {"NN", "E", "E", "2E", "E", "E", "2E", "E"};
    char tmp[2];
    char fits_err_msg[FITS_ERR_LEN];
    fitsfile *fptr;
    
    fits_open_table(&fptr, fname, READWRITE, &status);
    if (status) {
	fits_report_error(stderr, status);
	return NULL;
    }
    for (i = 0; i < 8; i++) {
	while ((colnum = col_exists(fptr, ttype[i]))) {
	    if (colnum < 0)
		return NULL;
	    fits_delete_col(fptr, colnum, &status);
	}
    }
    fits_get_num_rows(fptr, nrows, &status);
    fits_get_num_cols(fptr, &ncols, &status);
    if (status) {
       fits_close_file(fptr, &status);
       while (fits_read_errmsg(fits_err_msg))
	   msg_message(MSG_FATAL, 0, fits_err_msg);
       return NULL;
    }
    fits_get_colnum(fptr, 1, prefs.magname, &colnum, &status);
    fits_get_coltype(fptr, colnum, &type, &repeat, &width, &status);
    msg_message(MSG_DEBUG, 0, "%s in column %d has %ld elements", 
		prefs.magname, colnum, repeat);
    sprintf(tmp, "%ldE", repeat);
    tform[0] = strdup(tmp);

    fits_insert_cols(fptr, ncols + 1, 8, ttype, tform, &status);
    if (status) {
	msg_message(MSG_FATAL, 0, "Error while adding columns:");
        fits_close_file(fptr, &status);
	while (fits_read_errmsg(fits_err_msg))
	    msg_message(MSG_FATAL, 0, fits_err_msg);
       return NULL;
    }
    return fptr;
}


/*!\brief Read lensing catalog properties into arrays
 * \param fptr Open fits pointer to a table containing lensing properties
 * \param iter Number of iteration reading this catalog
 * \param blocksize Size of chunks to be read in each iteration
 * \param nrows Total number of rows in table
 * \param ind Array of galaxy indices used for matching
 * \param ora Array of observed right ascension
 * \param odec Array of observed declination
 * \param a00 Array of 0,0 component of lensing Jacobian
 * \param a01 Array of 0,1 component of lensing Jacobian
 * \param a10 Array of 1,0 component of lensing Jacobian
 * \param a11 Array of 1,1 component of lensing Jacobian
 * \return Number of items read, 0 in case of failure
 */
long
read_jacobian(fitsfile *fptr, long iter, long blocksize, long nrows, 
	      long *ind, double *ora, double *odec, double *a00, double *a01, 
	      double *a10, double *a11)
{
    long toread;
    char fits_err_msg[FITS_ERR_LEN];
    int status = 0;
    
    if (iter * blocksize > nrows)
	return 0;
    if ((iter + 1) * blocksize > nrows)
	toread = nrows - iter * blocksize;
    else
	toread = blocksize;

    fits_read_col(fptr, TLONG, 8, iter * blocksize + 1, 1, toread, NULL,
		  ind, NULL, &status);
    fits_read_col(fptr, TDOUBLE, 1, iter * blocksize + 1, 1, toread, NULL,
		  ora, NULL, &status);
    fits_read_col(fptr, TDOUBLE, 2, iter * blocksize + 1, 1, toread, NULL,
		  odec, NULL, &status);
    fits_read_col(fptr, TDOUBLE, 3, iter * blocksize + 1, 1, toread, NULL,
		  a00, NULL, &status);
    fits_read_col(fptr, TDOUBLE, 4, iter * blocksize + 1, 1, toread, NULL,
		  a01, NULL, &status);
    fits_read_col(fptr, TDOUBLE, 5, iter * blocksize + 1, 1, toread, NULL,
		  a10, NULL, &status);
    fits_read_col(fptr, TDOUBLE, 6, iter * blocksize + 1, 1, toread, NULL,
		  a11, NULL, &status);
    if (status) {
	msg_message(MSG_FATAL, 0, "Error reading shear catalog:");
	while (fits_read_errmsg(fits_err_msg))
	    msg_message(MSG_FATAL, 0, fits_err_msg);
	return 0;
    }
    return toread;
}


/*!\brief Determine whether an object has already been lensed, i.e., is a multiply imaged source
 * \param fptr Open fits file pointer to DES_mock table
 * \param idx Row number of the object to be checked
 * \return 1 if object was already lensed, 0 if not, -1 in case of error
 */
int
is_multiply_lensed(fitsfile *fptr, long idx)
{
    double val;
    int colnum;
    int status = 0;
    char fits_err_msg[FITS_ERR_LEN];

    /* We read the shear and surface mass density. If all  are zero, the
       object has not been lensed. */
    fits_get_colnum(fptr, 1, "KAPPA", &colnum, &status);
    fits_read_col(fptr, TDOUBLE, colnum, idx + 1, 1, 1, NULL, &val, NULL, 
		  &status);
    if (fabs(val) > 1e-5 && status == 0)
	return 1;
    fits_get_colnum(fptr, 1, "GAMMA", &colnum, &status);
    fits_read_col(fptr, TDOUBLE, colnum, idx + 1, 1, 1, NULL, &val, NULL, 
		  &status);
    if (fabs(val) > 1e-5 && status == 0)
	return 1;
    fits_read_col(fptr, TDOUBLE, colnum, idx + 1, 2, 1, NULL, &val, NULL, 
		  &status);
    if (fabs(val) > 1e-5 && status == 0)
	return 1;
    if (status) {
	msg_message(MSG_FATAL, 0, "Error writing Jacobian:");
        while (fits_read_errmsg(fits_err_msg))
            msg_message(MSG_FATAL, 0, fits_err_msg);
	fits_close_file(fptr, &status);
	return -1;
    }
    return 0;
}


/*!\brief Duplicate a multiply lensed object by copying the pertinent information to the end of the table
 * \param fptr Open fits file pointer to DES mock catalog table
 * \param nrows Number of rows the table currently has
 * \param idx Index position of the objects to be copied
 * \return The new number of rows, 0 in case of error
 */
long
duplicate_object(fitsfile *fptr, long nrows, long idx)
{
    long width;  /* row width in bytes */
    unsigned char *tmp;
    int status = 0;
    char fits_err_msg[FITS_ERR_LEN];    

    msg_message(MSG_INFO, 2, "Duplicating object %ld", idx);
    fits_read_key(fptr, TLONG, "NAXIS1", &width, NULL, &status);
    if (status) {
	msg_message(MSG_FATAL, 0, "Could not determine row width:");
	while (fits_read_errmsg(fits_err_msg))
	    msg_message(MSG_FATAL, 0, fits_err_msg);
	return 0;
    }
    CALLOC(tmp, width, "row copy buffer");
    fits_read_tblbytes(fptr, idx + 1, 1, width, tmp, &status);
    fits_write_tblbytes(fptr, nrows + 1, 1, width, tmp, &status);
    free(tmp);
    if (status) {
	msg_message(MSG_FATAL, 0, "Error copying rows:");
	while (fits_read_errmsg(fits_err_msg))
	    msg_message(MSG_FATAL, 0, fits_err_msg);
	return 0;
    }
    nrows++;
    fits_update_key(fptr, TLONG, "NAXIS2", &nrows, 0, &status);
    return nrows - 1;
}


/*!\brief Read lensed positions and Jacobian from lensing catalog and it lensing quantities (shear, surface mass density, magnification) to catalog
 * \param DES_mock Open fits file pointer to DES mock catalog table
 * \param fname File name of lensing catalog
 * \return Number of rows in case of success, 0 else
 */
long
add_lensing(fitsfile *DES_mock, char *fname)
{
    int colnum;
    long i, idx, iter, nread, nrows, blocksize, DES_nrows;
    long *ind;
    double *ora, *odec, *a00, *a01, *a10, *a11;
    double gamma[2], kappa, mu;
    int status = 0;
    char fits_err_msg[FITS_ERR_LEN];
    fitsfile *fptr;

    fits_get_num_rows(DES_mock, &DES_nrows, &status);
    fits_open_table(&fptr, fname, READONLY, &status);
    if (status) {
	fits_report_error(stderr, status);
	return 0;
    }
    fits_get_num_rows(fptr, &nrows, &status);
    fits_get_rowsize(fptr, &blocksize, &status);
    msg_message(MSG_DEBUG, 0, "Will read chunks of %ld rows", blocksize);
    CALLOC(ind, blocksize, "index");
    CALLOC(ora, blocksize, "observed RA");
    CALLOC(odec, blocksize, "observed DEC");
    CALLOC(a00, blocksize, "observed A00");
    CALLOC(a01, blocksize, "observed A01");
    CALLOC(a10, blocksize, "observed A10");
    CALLOC(a11, blocksize, "observed A11");
    iter = 0;
    while ((nread = read_jacobian(fptr, iter, blocksize, nrows, ind, ora, 
				  odec, a00, a01, a10, a11))) {
	msg_message(MSG_INFO, 1, "Reading chunk %ld/%ld", 
		    iter + 1, nrows / blocksize + 1);
	for (i = 0; i < nread; i++) {
	    idx = ind[i];
	    if (is_multiply_lensed(DES_mock, idx)) {
		idx = duplicate_object(DES_mock, DES_nrows, idx);
		if (idx == -1) {
		    msg_message(MSG_FATAL, 0,
			    "Error duplicating multiply lensed object.");
		    return 0;
		}
		DES_nrows = idx + 1;
	    }
	    gamma[0] = 0.5 * (a00[i] - a11[i]);
	    gamma[1] = -0.5 * (a01[i] + a10[i]);
	    kappa = 1 - 0.5 * (a00[i] + a11[i]);
	    mu = 1 / (a11[i] * a00[i] - a01[i] * a10[i]);

	    fits_get_colnum(DES_mock, CASEINSEN, "ORA", &colnum, &status);
	    fits_write_col(DES_mock, TDOUBLE, colnum, idx + 1, 1, 1, 
			   &ora[i], &status);
	    fits_get_colnum(DES_mock, CASEINSEN, "ODEC", &colnum, &status);
	    fits_write_col(DES_mock, TDOUBLE, colnum, idx + 1, 1, 1, 
			   &odec[i], &status);

	    fits_get_colnum(DES_mock, CASEINSEN, "GAMMA", &colnum, &status);
	    fits_write_col(DES_mock, TDOUBLE,  colnum, idx + 1,
			   1, 1, &gamma[0],  &status);
	    fits_write_col(DES_mock, TDOUBLE,  colnum, idx + 1,
			   2, 1, &gamma[1],  &status);
	    fits_get_colnum(DES_mock, CASEINSEN, "KAPPA", &colnum, &status);
	    fits_write_col(DES_mock, TDOUBLE, colnum, idx + 1,
			   1, 1, &kappa, &status);
	    fits_get_colnum(DES_mock, CASEINSEN, "MU", &colnum, &status);
	    fits_write_col(DES_mock, TDOUBLE, colnum, idx + 1,
			   1, 1, &mu, &status);

	    if (status) {
		msg_message(MSG_FATAL, 0, "Error writing Jacobian:");
		while (fits_read_errmsg(fits_err_msg))
		    msg_message(MSG_FATAL, 0, fits_err_msg);
		fits_close_file(DES_mock, &status);
		return 0;
	    }
	}
	iter++;
    }
    msg_message(MSG_DEBUG, 2, "Done");
    fits_close_file(fptr, &status);
    free(ind);
    free(ora);
    free(odec);
    free(a00);
    free(a01);
    free(a10);
    free(a11);
    return DES_nrows;
}




/*!\brief Read lensing quantities from catalog into arrays
 * \param fptr Fits file pointer to an open table
 * \param iter Number of iteration reading this catalog
 * \param blocksize Size of chunks to be read in each iteration
 * \param nrows Total number of rows in table
 * \param te Array of true (intrinsic) ellipticity
 * \param gamma Array of 2 component shear
 * \param kappa Aray of convergence
 * \return Number of items read, 0 in case of failure
 */
long
read_shear(fitsfile *fptr, long iter, long blocksize, long nrows, double *te,
	   double *gamma, double *kappa)
{
    long toread;
    int colnum;
    int status = 0;

    if (iter * blocksize > nrows)
	return 0;
    if ((iter + 1) * blocksize > nrows)
	toread = nrows - iter * blocksize;
    else
	toread = blocksize;

    fits_get_colnum(fptr, 1, prefs.ename, &colnum, &status);
    if (status) {
	fits_report_error(stderr, status);
	return 0;
    }
    fits_read_col(fptr, TDOUBLE, colnum, iter * blocksize + 1, 1, 2 * toread, 
		  NULL, te, NULL, &status);
    fits_get_colnum(fptr, 1, "GAMMA", &colnum, &status);
    if (status) {
	fits_report_error(stderr, status);
	return 0;
    }
    fits_read_col(fptr, TDOUBLE, colnum, iter * blocksize + 1, 1, 2 * toread, 
		  NULL, gamma, NULL, &status);
    fits_get_colnum(fptr, 1, "KAPPA", &colnum, &status);
    if (status) {
	fits_report_error(stderr, status);
	return 0;
    }
    fits_read_col(fptr, TDOUBLE, colnum, iter * blocksize + 1, 1, toread, NULL,
		  kappa, NULL, &status);
    return toread;
}


/*!\brief Apply magnification to object sizes
 *\param fptr pointer to a fitsfile with size and magnification column
 *\returns 1 in case of success, 0 in case of error
 */
int
magnify_size(fitsfile *fptr)
{
    int colnum;
    char fits_err_msg[FITS_ERR_LEN];
    char calcstr[1024];
    int status = 0;

    fits_get_colnum(fptr, 1, prefs.sname, &colnum, &status);
    if (status == COL_NOT_FOUND) {
	    msg_message(MSG_NOTICE, 0, 
			"Column %s not found. Is this a magnitude only file?", 
			prefs.sname);
	    msg_message(MSG_WARNING, 0, 
			"Will proceed but will not attempt to lens sizes/shapes");
	    prefs.magonly = 1;
	    status = 0;
    } else {
	msg_message(MSG_INFO, 0, "Applying magnification to size ...");
	snprintf(calcstr, 1024, "%s * sqrt(abs(MU))", prefs.sname);
	fits_calculator(fptr, calcstr, fptr, prefs.lenssname, NULL, &status);
	if (status) {
	    fits_close_file(fptr, &status);
	    msg_message(MSG_FATAL, 0, "Error computing %:", prefs.lenssname);
	    while (fits_read_errmsg(fits_err_msg))
		msg_message(MSG_FATAL, 0, fits_err_msg);
	    return 0;
	}
	msg_message(MSG_INFO, 0, "Done");
    }
    return 1;
}


/*!\brief Apply magnification to magnitude column
 *\param fptr pointer to a fitsfile with magnitude and magnification column
 *\returns 1 in case of success, 0 in case of error
 */
int
magnify_magnitudes(fitsfile *fptr)
{
    char fits_err_msg[FITS_ERR_LEN];
    char calcstr[1024];
    int status = 0;

    msg_message(MSG_INFO, 0, "Applying magnification to magnitudes ...");
    snprintf(calcstr, 1024, "%s - 2.5 *  log10(abs(MU))", prefs.magname);
    fits_calculator(fptr, calcstr, fptr, prefs.lensmagname, NULL, &status);
    if (status) {
        fits_close_file(fptr, &status);
        msg_message(MSG_FATAL, 0, "Error computing %s:", prefs.lensmagname);
        while (fits_read_errmsg(fits_err_msg))
            msg_message(MSG_FATAL, 0, fits_err_msg);
        return 0;
    }
    msg_message(MSG_INFO, 0, "Done.");
    return 1;
}


/*!\brief Apply reduced shear to mock shapes of catalogs
 *\param fptr pointer to fits file with intrinsic shape and shear/kappa
 *\param nrows Number of rows in input table
 *\returns 1 in case of success, 0 in case of error
 */
int 
lens_shapes(fitsfile *fptr, long nrows)
{
    double *te, *gamma, *kappa;
    gsl_complex epss, eps, g;
    int colnum;
    long i, nread, iter, blocksize;
    char fits_err_msg[FITS_ERR_LEN];
    int status = 0;

    fits_get_rowsize(fptr, &blocksize, &status);
    fits_get_colnum(fptr, 1, prefs.shearname, &colnum, &status);
    if (status) {
	fits_report_error(stderr, status);
	return 0;
    }
    CALLOC(te, 2 * blocksize, "te");
    CALLOC(gamma, 2 * blocksize, "gamma");
    CALLOC(kappa, blocksize, "kappa");
    iter = 0;
    msg_message(MSG_INFO, 0, "Shearing galaxies ...");
    while ((nread = read_shear(fptr, iter, blocksize, nrows, te, gamma, 
			       kappa))) {
	msg_message(MSG_INFO, 1, "Lensing chunk %ld/%ld", 
		    iter + 1, nrows / blocksize + 1);
	for (i = 0; i < nread; i++) {
	    epss = gsl_complex_rect(te[2 * i], te[2 * i+1]);
	    g = gsl_complex_div(gsl_complex_rect(gamma[2 * i], 
						 gamma[2 * i + 1]),
				gsl_complex_rect(1 - kappa[i], 0));
	    if (gsl_complex_abs(g) <= 1) {
		eps = gsl_complex_div(gsl_complex_add(epss, g), 
				      gsl_complex_add(GSL_COMPLEX_ONE, 
						      gsl_complex_mul(epss,
								      gsl_complex_conjugate(g))));
	    } else {
		eps = gsl_complex_div(gsl_complex_add(GSL_COMPLEX_ONE, 
						      gsl_complex_mul(g, 
								      gsl_complex_conjugate(epss))),
				      gsl_complex_add(
						      gsl_complex_conjugate(epss),
						      gsl_complex_conjugate(g)));
	    }
	    if (gsl_complex_abs(eps) > 1 || gsl_complex_abs(epss) > 1) {
		msg_message(MSG_FATAL, 0, 
			    "|eps| = %.10g > 1 OR |eps_s| = %.10g > 1", 
			    gsl_complex_abs(eps), gsl_complex_abs(epss));
		msg_message(MSG_INFO, 0, "eps_s = %.10g + %.10gj", 
			    GSL_REAL(epss), GSL_IMAG(epss));
		msg_message(MSG_INFO, 0, "g     = %.10g + %.10gj", 
			    GSL_REAL(g), GSL_IMAG(g));
		msg_message(MSG_INFO, 0, "eps   = %.10g + %.10gj", 
			    eps.dat[0], eps.dat[1]);
		return 0;
	    }
	    
	    fits_write_col(fptr, TDOUBLE, colnum, iter * blocksize + i + 1,
			   1, 1, &eps.dat[0], &status);
	    fits_write_col(fptr, TDOUBLE, colnum, iter * blocksize + i + 1,
			   2, 1, &eps.dat[1], &status);
	    if (status) {
		while (fits_read_errmsg(fits_err_msg))
		    msg_message(MSG_FATAL, 0, fits_err_msg);
		return 0;
	    }
	}
	iter++;
    }
    if (iter < nrows / blocksize) {
	/* read_shear encountered an error before we reached the end */
	fits_report_error(stderr, status);
	return 0;
    }
    msg_message(MSG_INFO, 2, "Done.");
    free(te);
    free(gamma);
    free(kappa);
    return 1;
}


/*!\brief Apply lensing Jacobian to objects, magnify magnitude and size, shear galaxies
 * \param fptr fitsfile pointer to an open table
 * \param nrows number of rows in fits table
 * \returns 1 in case of success, 0 else
 */
int
lens_galaxies(fitsfile *fptr, long nrows)
{
    if (!prefs.magonly)
	if (!magnify_size(fptr))
	    return 0;
    
    if (!magnify_magnitudes(fptr))
	return 0;

    if (!prefs.magonly)
	if (!lens_shapes(fptr, nrows))
	    return 0;

    return 1;
}

/*!\brief Delete all lensing columns that are not magnitude information. This is called for non-shape catalogs.
 *\params fptr Pointer to output fits table
 */
void
delete_nonmaginfo(fitsfile *fptr)
{
    int i;
    int colnum;
    char *ttype[] = {"ORA", "ODEC", "GAMMA", "KAPPA", "MU", 
		     prefs.shearname, prefs.lenssname};
    int status = 0;

    for (i = 0; i < 7; i++) {
	fits_get_colnum(fptr, CASEINSEN, ttype[i], &colnum, &status);
	if (!status) {
	    fits_delete_col(fptr, colnum, &status);
	    msg_message(MSG_INFO, 0, "Deleted column %s", ttype[i]);
	}
    }
}


int 
main(int argc, char **argv)
{
    long nrows;
    int status = 0;
    unsigned short havegal = 0;
    fitsfile *DES_mock;
    int a, narg, opt;
    static char prefsname[PATH_MAX];
    char **argkey, **argval;

    msg_init(MSG_INFO, NULL);
    if (argc < 2) {
        msg_message(MSG_INFO, 0, "%s-%s", "apply-lensing", VERSION);
        msg_message(MSG_INFO, 0, "Written by %s", COPYRIGHT);
        msg_message(MSG_FATAL, 0, "%s", "SYNTAX:");
        msg_message(MSG_FATAL, 0, 
		    "    apply-lensing DES-mock.fits shearcat.fits");
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
                default: {
                    msg_message(MSG_FATAL, 0, "SYNTAX:\n%s", SYNTAX);
		    msg_message(MSG_FATAL, 0, 
				"apply-lensing DES-mock.fits shearcat.fits");
		}
                }
            } else {
                argkey[narg] = &argv[a][1];
                argval[narg++] = argv[++a];
            }
        } else {
	    if (!havegal) {
		snprintf(prefs.incat, PATH_MAX, "%s", argv[a]);
		havegal = 1;
	    }
	    else
		snprintf(prefs.shearcat, PATH_MAX, "%s", argv[a]);
        }
    }

    readprefs(prefsname, argkey, argval, narg);

    msg_init(prefs.verbose, NULL);
    msg_message(MSG_INFO, 0, "%s-%s", "apply-lensing", VERSION);
    msg_message(MSG_INFO, 0, "by %s", COPYRIGHT);
    msg_message(MSG_DEBUG, 0, "Verbosity level: %d", prefs.verbose);

    msg_message(MSG_INFO, 0, "Adding output columns to %s", prefs.incat);
    msg_message(MSG_INFO, 0, "Reading shear from %s", prefs.shearcat);
    if (!(DES_mock = prepare_table(prefs.incat, &nrows)))
	exit(EXIT_FAILURE);
    msg_message(MSG_INFO, 0, "Done.");
    fits_flush_file(DES_mock, &status);
    msg_message(MSG_INFO, 0, 
		"Will add lensed positions and lensing quantities ...");
    if (!(nrows = add_lensing(DES_mock, prefs.shearcat)))
	exit(EXIT_FAILURE);
    msg_message(MSG_INFO, 0, "Done.");
    fits_flush_file(DES_mock, &status);
    msg_message(MSG_INFO, 0, "Will apply lensing to galaxies ...");
    if (!lens_galaxies(DES_mock, nrows))
	exit(EXIT_FAILURE);
    msg_message(MSG_INFO, 0, "Done.");
    if (prefs.magonly)
	delete_nonmaginfo(DES_mock);
    fits_close_file(DES_mock, &status);
    if (status) {
	fits_report_error(stderr, status);
	exit(EXIT_FAILURE);
    }

    msg_finish();
    return EXIT_SUCCESS;
}
