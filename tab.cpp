#include <bits/stdc++.h>

using namespace std;

int main() {
	for (int i = 1; i <= 1000000; i++) {
		cout << unsigned(-1) / i + 1 << ",\n"[i == 1000000];
	}
}