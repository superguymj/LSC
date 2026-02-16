#include <vector>

template <class T>
struct Array2D {
	int n;
	std::vector<T> a;
	Array2D() {}
	Array2D(int n) : n(n), a(n * n) {}

	T& operator()(int i, int j) {
		return a[i * n + j];
	}

	const T& operator()(int i, int j) const {
		return a[i * n + j];
	}
};
template <class T>
struct Array3D {
	int n;
	std::vector<T> a;
	Array3D() {}
	Array3D(int n) : n(n), a(n * n * n) {}

	T& operator()(int i, int j, int k) {
		return a[i * n * n + j * n + k];
	}
	const T& operator()(int i, int j, int k) const {
		return a[i * n * n + j * n + k];
	}
};