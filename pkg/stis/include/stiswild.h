#ifndef INCL_STISWILD_H
#define INCL_STISWILD_H

/* This file defines the wildcards for use when comparing values
   (e.g. from the input header) with values read from a table row.
   The IGNORE values mean that the column is not relevant for the
   particular table.
*/

# define  INT_WILDCARD      -1
# define  STRING_WILDCARD  "ANY"

# define  INT_IGNORE        -999
# define  STRING_IGNORE_1  "N/A"
# define  STRING_IGNORE_2  "NOT APPLICABLE"

#endif /* INCL_STISWILD_H */
