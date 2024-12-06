AM_CPPFLAGS = -iquote $(top_builddir) -iquote $(srcdir)/.. -iquote $(top_srcdir) -iquote $(top_srcdir)/cwds

noinst_LTLIBRARIES = libmath.la

libmath_la_SOURCES = \
	AnalyzedCubic.cxx \
	bracket_zero.cxx \
	CubicPolynomial.cxx \
	Direction.cxx \
	Line.cxx \
	LinePiece.cxx \
	Point.cxx \
	Polynomial.cxx \
	AnalyzedCubic.h \
	bracket_zero.h \
	CubicPolynomial.h \
	Direction.h \
	Line.h \
	LinePiece.h \
	Point.h \
	Polynomial.h \
	QuadraticPolynomial.h \
	Vector.h

libmath_la_CXXFLAGS = @LIBCWD_R_FLAGS@

# --------------- Maintainer's Section

if MAINTAINER_MODE
distclean-local:
	rm -f *.s *.ii

clean-local:
endif

MAINTAINERCLEANFILES = $(srcdir)/Makefile.in
