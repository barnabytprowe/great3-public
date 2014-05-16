inputname$ = $PATHNAME() + "\GREAT3_COEFFICIENTS\ZEMAXInput.dat"
OPEN inputname$
config$ = $PATHNAME() + "\des-proto.CFG"
PRINT config$
LABEL 1
READ x, y
IF (EOFF()) THEN GOTO 2
PRINT x, " ", y
SETSYSTEMPROPERTY 102, 9, x
SETSYSTEMPROPERTY 103, 9, y
outputname$ = $PATHNAME() + "\GREAT3_COEFFICIENTS\" +  $STR(x) + "_" + $STR(y) + ".txt"
GETTEXTFILE outputname$, Zst, config$, 1
PRINT outputname$
GOTO 1
LABEL 2
CLOSE