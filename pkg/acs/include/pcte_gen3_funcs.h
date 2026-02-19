/* These functions are for dopcte-gen3.c in ACSCTE.
   They are not meant to be used outside that module and they were static.
   We are only exposing them here so they can be tested in CI. */

int rotateAmpData_acscte(FloatTwoDArray * amp, const unsigned ampID);
int derotateAmpData_acscte(FloatTwoDArray * amp, const unsigned ampID);
