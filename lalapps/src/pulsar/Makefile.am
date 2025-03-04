SUBDIRS = \
	Tools \
	SFTTools \
	FITSTools \
	MakeData \
	Fscan \
	CreateEphemeris \
	Fstatistic \
	HeterodyneSearch \
	Hough \
	HoughFstat \
	GCT \
	CrossCorr \
	SidebandSearch \
	TwoSpect \
	Weave \
	$(END_OF_LINE)

# Because many tests in lalapps/src/pulsar/ call executables from other
# subdirectories in lalapps/src/pulsar/, it is safest to make sure the
# whole of src/pulsar is built first
check-recursive: all-recursive

# All programs under lalapps/src/pulsar/ [1] must be executed by at least one
# test script run by "make check". This applies to C executables, Python
# scripts, etc. The following rule enforces this requirement.
#
# [1] Certain directories are currently exempt from this requirement; in these
# directories, Makefile.am includes lalapps.am instead of lalapps_pulsar.am
.PHONY: check-execs-are-tested
check-local: check-execs-are-tested
check-execs-are-tested:
	@tested_logfiles=`find $(builddir) -name '.tested_*.log'`; \
	not_tested=`$(AWK) ' \
	$$1 == "programs" { for (i = 3; i <= NF; ++i) { prog[$$i] = 1; progpath[$$i] = $$2 "/" $$i } } \
	$$1 == "tested" { for (i = 2; i <= NF; ++i) { tested[$$i] = 1; } } \
	END { if (prog) { for (e in prog) { if (!tested[e]) { print progpath[e]; } } } else { print "NONE"; } } \
	' $${tested_logfiles} /dev/null | $(SED) 's|^$(subdir)/||' | sort`; \
	if test "X$${not_tested}" = XNONE; then \
		exit 0; \
	fi; \
	if test "X$${not_tested}" != X; then \
		printf "\nERROR: the following executables under lalapps_pulsar.am test harness are not tested:\n" >&2; \
		printf "   %s\n" $${not_tested} >&2; \
		printf "\n" >&2; \
		exit 1; \
	fi; \
	printf "OK: all executables under lalapps_pulsar.am test harness are tested\n"
