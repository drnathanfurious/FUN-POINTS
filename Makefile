all:
	ifort -m64 -r8 -i8 -fast -opt_report -o out main.f90
