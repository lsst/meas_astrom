#include "starkd.h"

static PyObject* pyval_int(void* v, int index) {
	PyObject* pyval;
	int* i = (int*)v;
	pyval = PyInt_FromLong((long)i[index]);
	assert(pyval);
	return pyval;
}
static PyObject* pyval_int64(void* v, int index) {
	PyObject* pyval;
	int64_t* i = (int64_t*)v;
	pyval = PyLong_FromLongLong((long long)i[index]);
	assert(pyval);
	return pyval;
}
static PyObject* pyval_double(void* v, int index) {
	PyObject* pyval;
	double* i = (double*)v;
	pyval = PyFloat_FromDouble(i[index]);
	assert(pyval);
	return pyval;
}

static PyObject* array_to_pylist(void* X, int N, PyObject* (*pyval)(void* v, int i)) {
	PyObject* pylist;
	int i;
	pylist = PyList_New(N);
	assert(pylist);
	for (i=0; i<N; i++) {
		PyObject* v;
		v = pyval(X, i);
		assert(v);
		PyList_SetItem(pylist, i, v);
	}
	return pylist;
}
static PyObject* array_to_pylist2(void* vdata, int N, int D, PyObject* (*pyval)(void* v, int i)) {
	PyObject* pylist;
	int i,j;
	pylist = PyList_New(N);
	assert(pylist);
	for (i=0; i<N; i++) {
		PyObject* v;
		PyObject* row = PyList_New(D);
		assert(row);
		for (j=0; j<D; j++) {
			v = pyval(vdata, i*D+j);
			assert(v);
			PyList_SetItem(row, j, v);
		}
		PyList_SetItem(pylist, i, row);
	}
	return pylist;
}

static PyObject* arrayd_to_pylist2(double* xyz, int N, int D) {
	return array_to_pylist2(xyz, N, D, pyval_double);
}
static PyObject* arrayi_to_pylist(int* X, int N) {
	return array_to_pylist(X, N, pyval_int);
}

/*
//static 
PyObject* starkd_search_stars(PyObject* self, PyObject* args) {
	startree_t* s;
	double ra, dec, radius;
    int N;
	PyObject* pyxyz;
	PyObject* pyradec;
	PyObject* pyinds;
	double* xyzres = NULL;
	double* radecres = NULL;
	int* inds = NULL;
	unsigned char tag;
	int i, C;
	PyObject* pydict;
    if (!PyArg_ParseTuple(args, "ldddb", &s, &ra, &dec, &radius, &tag)) {
        PyErr_SetString(PyExc_ValueError, "need four args: starkd, ra, dec, radius");
        return NULL;
	}
 */

static PyObject* starkd_return_pystuff(startree_t* s, double* radecres, double* xyzres, int* inds, int N) {
	unsigned char tag = 1;
	PyObject* pyxyz;
	PyObject* pyradec;
	PyObject* pyinds;
	int i, C;
	PyObject* pydict;

	pyxyz = arrayd_to_pylist2(xyzres, N, 3);
	pyradec = arrayd_to_pylist2(radecres, N, 2);
	pyinds = arrayi_to_pylist(inds, N);
	if (!tag)
		return Py_BuildValue("(OOO)", pyxyz, pyradec, pyinds);
	if (!startree_has_tagalong(s) || N == 0)
		return Py_BuildValue("(OOOO)", pyxyz, pyradec, pyinds, PyDict_New());
	C = startree_get_tagalong_N_columns(s);
	pydict = PyDict_New();
	for (i=0; i<C; i++) {
		PyObject* pyval;
		void* vdata;
		const char* name;
		tfits_type ft, readtype;
		int arr;
		PyObject* (*pyvalfunc)(void* v, int i);
		name = startree_get_tagalong_column_name(s, i);
		//ft = startree_get_tagalong_column_fits_type(s, i);
		// v0.30:
		//ft = fitstable_get_fits_column_type(startree_get_tagalong(s), i);
		// HACK...
		ft = startree_get_tagalong(s)->table->col[i].atom_type;

		readtype = ft;
		switch (ft) {
			/*
			 // FIXME -- special-case this?
			 case TFITS_BIN_TYPE_A:
			 case TFITS_BIN_TYPE_B:
			 case TFITS_BIN_TYPE_X:
			 // these are actually the characters 'T' and 'F'.
			 case TFITS_BIN_TYPE_L:
			 pytype = NPY_BYTE;
			 break;
			 */
		case TFITS_BIN_TYPE_D:
		case TFITS_BIN_TYPE_E:
			readtype = TFITS_BIN_TYPE_D;
			pyvalfunc = pyval_double;
			break;
		case TFITS_BIN_TYPE_I:
		case TFITS_BIN_TYPE_J:
			readtype = TFITS_BIN_TYPE_J;
			pyvalfunc = pyval_int;
			break;
		case TFITS_BIN_TYPE_K:
			readtype = TFITS_BIN_TYPE_K;
			pyvalfunc = pyval_int64;
			break;
		default:
			PyErr_Format(PyExc_ValueError, "failed to map FITS type %i to numpy type, for column \"%s\"", ft, name);
			return NULL;
		}
		vdata = fitstable_read_column_array_inds(startree_get_tagalong(s), name, readtype, inds, N, &arr);
		if (!vdata) {
			PyErr_Format(PyExc_ValueError, "fail to read tag-along column \"%s\"", name);
			return NULL;
		}
		if (arr > 1)
			pyval = array_to_pylist2(vdata, N, arr, pyvalfunc);
		else
			pyval = array_to_pylist(vdata, N, pyvalfunc);
		if (PyDict_SetItemString(pydict, name, pyval)) {
			PyErr_Format(PyExc_ValueError, "fail to set tag-along column value, for \"%s\"", name);
			return NULL;
		}
	}
	return Py_BuildValue("(OOOO)", pyxyz, pyradec, pyinds, pydict);
}

PyObject* starkd_search_stars_in_field(startree_t* s, tan_t* tanwcs, double pixelmargin) {
    int N;
	double* xyzres = NULL;
	double* radecres = NULL;
	int* inds = NULL;
	int i;
	int Nkeep;

	{
		double xyz[3];
		double r2;
		double px = (tanwcs->imagew + 1.0) / 2.0;
		double py = (tanwcs->imageh + 1.0) / 2.0;
		tan_pixelxy2xyzarr(tanwcs, px, py, xyz);
		r2 = arcsec2distsq((hypot(px, py) + pixelmargin) * tan_pixel_scale(tanwcs));
		startree_search_for(s, xyz, r2, &xyzres, &radecres, &inds, &N);
		//printf("Found %i index stars\n", N);
	}

	assert(N == 0 || xyzres);
	assert(N == 0 || radecres);
	assert(N == 0 || inds);

	Nkeep = 0;
	for (i=0; i<N; i++) {
		double x,y;
		if (!tan_radec2pixelxy(tanwcs, radecres[2*i+0], radecres[2*i+1], &x, &y))
			continue;
		if (x < pixelmargin || y < pixelmargin ||
			x > tanwcs->imagew + pixelmargin || y > tanwcs->imageh + pixelmargin)
			continue;

		if (Nkeep != i) {
			memcpy(radecres + 2*Nkeep, radecres + 2*i, 2*sizeof(double));
			memcpy(xyzres   + 3*Nkeep, xyzres   + 3*i, 3*sizeof(double));
			inds[Nkeep] = inds[i];
		}
		Nkeep++;
	}
	// FIXME -- realloc?

	return starkd_return_pystuff(s, radecres, xyzres, inds, Nkeep);
}

PyObject* starkd_search_stars(startree_t* s, double ra, double dec, double radius) {
    int N;
	double* xyzres = NULL;
	double* radecres = NULL;
	int* inds = NULL;

	//startree_search_for_radec(s, ra, dec, radius, &xyzres, &radecres, &inds, &N);
	// v0.30:
	{
		double xyz[3];
		double r2;
		radecdeg2xyzarr(ra, dec, xyz);
		r2 = deg2distsq(radius);
		startree_search_for(s, xyz, r2, &xyzres, &radecres, &inds, &N);
		//printf("Found %i index stars\n", N);
	}

	assert(N == 0 || xyzres);
	assert(N == 0 || radecres);
	assert(N == 0 || inds);

	return starkd_return_pystuff(s, radecres, xyzres, inds, N);
}
