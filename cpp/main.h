/* Code generated by cmd/cgo; DO NOT EDIT. */

/* package command-line-arguments */


#line 1 "cgo-builtin-export-prolog"

#include <stddef.h>

#ifndef GO_CGO_EXPORT_PROLOGUE_H
#define GO_CGO_EXPORT_PROLOGUE_H

#ifndef GO_CGO_GOSTRING_TYPEDEF
typedef struct { const char *p; ptrdiff_t n; } _GoString_;
#endif

#endif

/* Start of preamble from import "C" comments.  */




/* End of preamble from import "C" comments.  */


/* Start of boilerplate cgo prologue.  */
#line 1 "cgo-gcc-export-header-prolog"

#ifndef GO_CGO_PROLOGUE_H
#define GO_CGO_PROLOGUE_H

typedef signed char GoInt8;
typedef unsigned char GoUint8;
typedef short GoInt16;
typedef unsigned short GoUint16;
typedef int GoInt32;
typedef unsigned int GoUint32;
typedef long long GoInt64;
typedef unsigned long long GoUint64;
typedef GoInt64 GoInt;
typedef GoUint64 GoUint;
typedef size_t GoUintptr;
typedef float GoFloat32;
typedef double GoFloat64;
#ifdef _MSC_VER
#include <complex.h>
typedef _Fcomplex GoComplex64;
typedef _Dcomplex GoComplex128;
#else
typedef float _Complex GoComplex64;
typedef double _Complex GoComplex128;
#endif

/*
  static assertion to make sure the file is being used on architecture
  at least with matching size of GoInt.
*/
typedef char _check_for_64_bit_pointer_matching_GoInt[sizeof(void*)==64/8 ? 1:-1];

#ifndef GO_CGO_GOSTRING_TYPEDEF
typedef _GoString_ GoString;
#endif
typedef void *GoMap;
typedef void *GoChan;
typedef struct { void *t; void *v; } GoInterface;
typedef struct { void *data; GoInt len; GoInt cap; } GoSlice;

#endif

/* End of boilerplate cgo prologue.  */

#ifdef __cplusplus
extern "C" {
#endif


/* Return type for ENU2AER */
struct ENU2AER_return {
	GoFloat64 az; /* az */
	GoFloat64 el; /* el */
	GoFloat64 srange; /* srange */
};
extern __declspec(dllexport) struct ENU2AER_return ENU2AER(GoFloat64 e, GoFloat64 n, GoFloat64 u);

/* Return type for AER2ENU */
struct AER2ENU_return {
	GoFloat64 e; /* e */
	GoFloat64 n; /* n */
	GoFloat64 u; /* u */
};
extern __declspec(dllexport) struct AER2ENU_return AER2ENU(GoFloat64 az, GoFloat64 el, GoFloat64 srange);

/* Return type for AER2ECEF */
struct AER2ECEF_return {
	GoFloat64 x; /* x */
	GoFloat64 y; /* y */
	GoFloat64 z; /* z */
};
extern __declspec(dllexport) struct AER2ECEF_return AER2ECEF(GoFloat64 az, GoFloat64 el, GoFloat64 srange, GoFloat64 lat0, GoFloat64 lon0, GoFloat64 h0, GoString coordinateSystem);

/* Return type for ECEF2AER */
struct ECEF2AER_return {
	GoFloat64 az; /* az */
	GoFloat64 el; /* el */
	GoFloat64 srange; /* srange */
};
extern __declspec(dllexport) struct ECEF2AER_return ECEF2AER(GoFloat64 x, GoFloat64 y, GoFloat64 z, GoFloat64 lat0, GoFloat64 lon0, GoFloat64 h0, GoString coordinateSystem);

/* Return type for Geodetic2ECEF */
struct Geodetic2ECEF_return {
	GoFloat64 x; /* x */
	GoFloat64 y; /* y */
	GoFloat64 z; /* z */
};
extern __declspec(dllexport) struct Geodetic2ECEF_return Geodetic2ECEF(GoFloat64 lat, GoFloat64 lon, GoFloat64 alt, GoString coordinateSystem);

/* Return type for ECEF2Geodetic */
struct ECEF2Geodetic_return {
	GoFloat64 lat; /* lat */
	GoFloat64 lon; /* lon */
	GoFloat64 alt; /* alt */
};
extern __declspec(dllexport) struct ECEF2Geodetic_return ECEF2Geodetic(GoFloat64 x, GoFloat64 y, GoFloat64 z, GoString coordinateSystem);

/* Return type for ECEF2ENU */
struct ECEF2ENU_return {
	GoFloat64 e; /* e */
	GoFloat64 n; /* n */
	GoFloat64 u; /* u */
};
extern __declspec(dllexport) struct ECEF2ENU_return ECEF2ENU(GoFloat64 x, GoFloat64 y, GoFloat64 z, GoFloat64 lat0, GoFloat64 lon0, GoFloat64 h0, GoString coordinateSystem);

/* Return type for ENU2ECEF */
struct ENU2ECEF_return {
	GoFloat64 x; /* x */
	GoFloat64 y; /* y */
	GoFloat64 z; /* z */
};
extern __declspec(dllexport) struct ENU2ECEF_return ENU2ECEF(GoFloat64 e, GoFloat64 n, GoFloat64 u, GoFloat64 lat0, GoFloat64 lon0, GoFloat64 h0, GoString coordinateSystem);
extern __declspec(dllexport) GoFloat64 juliandate(GoInt64 timestamp);
extern __declspec(dllexport) GoFloat64 greenwichsrt(GoFloat64 jd);

/* Return type for ECI2ECEF */
struct ECI2ECEF_return {
	GoFloat64 xEcef; /* xEcef */
	GoFloat64 yEcef; /* yEcef */
	GoFloat64 zEcef; /* zEcef */
};
extern __declspec(dllexport) struct ECI2ECEF_return ECI2ECEF(GoFloat64 x, GoFloat64 y, GoFloat64 z, GoInt64 timestamp);

/* Return type for ECEF2ECI */
struct ECEF2ECI_return {
	GoFloat64 xEci; /* xEci */
	GoFloat64 yEci; /* yEci */
	GoFloat64 zEci; /* zEci */
};
extern __declspec(dllexport) struct ECEF2ECI_return ECEF2ECI(GoFloat64 x, GoFloat64 y, GoFloat64 z, GoInt64 timestamp);

#ifdef __cplusplus
}
#endif
