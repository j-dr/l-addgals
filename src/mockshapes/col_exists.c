#include <fitsio.h>

#include "define.h"
#include "message.h"

/*!\brief Check whether a column in a table already exists
 * \param fptr fitsfile pointer to an open table
 * \param colname name of the column whose existence is to be checked
 * \return Column number if column exists, -1 in case of error, 0 otherwise
 */
int
col_exists(fitsfile *fptr, char *colname)
{
    int colnum;
    char fits_err_msg[FITS_ERR_LEN];
    int status = 0;

    fits_get_colnum(fptr, CASEINSEN, colname, &colnum, &status);
    if (status == 0 || status == COL_NOT_UNIQUE) {
	msg_message(MSG_WARNING, 0, 
		    "Output column %s already exists. Will be deleted.",
		    colname);
	return colnum;
    } else if (status != COL_NOT_FOUND) {
	fits_close_file(fptr, &status);
	msg_message(MSG_FATAL, 0, 
		    "Could not determine whether column %s already exists.", 
		    colname);
        while (fits_read_errmsg(fits_err_msg))
            msg_message(MSG_FATAL, 0, fits_err_msg);
	return -1;
    }
    /* Empty the CFITSIO error stack from the above calls */
    while (fits_read_errmsg(fits_err_msg)) 
	msg_message(MSG_DEBUG, 0, fits_err_msg);
    
    return 0;
}
