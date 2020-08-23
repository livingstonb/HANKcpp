#ifndef _HANK_BOOST_H
#define _HANK_BOOST_H

#include <boost/multi_array.hpp>

template <typename T, size_t N>
using boost_array_type = boost::multi_array<T, N>;

template <typename T, size_t N>
using boost_array_shape = boost::array<typename boost_array_type<T, N>::index, N>;

using boost3d = boost_array_type<double, 3>;

using boost1d = boost_array_type<double, 1>;

using boost3dshape = boost_array_shape<double, 3>;

using boost1dshape = boost_array_shape<double, 1>;

using range = boost::multi_array_types::index_range;

using boost_index = boost3d::index_gen;

template<typename T, size_t N>
boost_array_type<T, N> new_array(const boost_array_shape<T, N>& shape) {
	boost_array_type<T, N> arr(shape);

	return arr;
}

template<typename T, size_t N>
boost_array_type<T, N> reshape_array(const boost_array_type<T, N>& arr, const boost_array_shape<T, N>& shape) {
	boost_array_type<T, N> out = arr;
	out.reshape(shape);
	return out;
}

template<typename T>
boost_array_type<T, 3> flatten_array3d(const boost_array_type<T, 3>& arr) {
	boost_array_type<T, 3> out = reshape_array(arr, {static_cast<int>(arr.num_elements()), 1, 1});
	return out;
}


#endif