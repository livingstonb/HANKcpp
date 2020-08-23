
uint subdim0;
uint subdim1;
uint subdim2;
void set_dims_3d(std::vector<uint> dims) {
	if ( dims.size() >= 1 )
		subdim0 = dims[0];

	if ( dims.size() >= 2 )
		subdim1 = dims[1];

	if ( dims.size() >= 3 )
		subdim2 = dims[2];
}
inline Scalar as3d(uint i, uint j, uint k) const {
	return this->operator()(i + subdim0 * j, k);
}

inline Scalar& as3d(uint i, uint j, uint k) {
	return this->operator()(i + subdim0 * j, k);
}