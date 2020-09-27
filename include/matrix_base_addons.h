
uint subdim0;
uint subdim1;
uint subdim2;
void set_dims_3d(const int* dims, int n) {
	if ( n >= 1 )
		subdim0 = dims[0];

	if ( n >= 2 )
		subdim1 = dims[1];

	if ( n >= 3 )
		subdim2 = dims[2];
}

inline Scalar as3d(uint i, uint j, uint k) const {
	return this->operator()(i + subdim0 * j, k);
}

inline Scalar& as3d(uint i, uint j, uint k) {
	return this->operator()(i + subdim0 * j, k);
}

// template<typename ScalarThis, typename ScalarR>
// struct auto_cast


// operator