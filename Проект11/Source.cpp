#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <vector>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <map>
#include <ctime>
#define pi 3.14159265358979323
using namespace std;

#define kmp_out

#define is_honeycomb
#define kohon_map_out
double dist(const vector<double>&a, const vector<double>&b, bool is_sqrt) {
	double res = 0;
	for (int i = 0; i < a.size(); ++i) {
		res += (a[i] - b[i]) * (a[i] - b[i]);
	}
	if (is_sqrt)
		res = sqrt(res);
	return res;
}

struct kohonen_map {
	int dim, width, height, cnt;
	vector<vector<double> > w;
	vector<vector<double>> x;
	kohonen_map(int _dim, int _w, int _h, vector<pair<double, double> >& border) {
		dim = _dim;
		width = _w;
		height = _h;
		//init honeycomb
		for (int i = 0; i < height; ++i) {
			double d = 0;
#ifdef is_honeycomb
			d = (i % 2) * 0.5;
#endif
			for (int j = 0; j < width; ++j) {
				x.push_back({ (double)i, (double)j + d });
			}
		}
		//init weights
		cnt = _w * _h;
		w.resize(cnt, vector<double>(dim, 0));
		for (int i = 0; i < cnt; ++i) {
			for (int j = 0; j < dim; ++j) {
				w[i][j] = ((double)rand() / RAND_MAX) * (border[j].second - border[j].first) + border[j].first;
			}
		}
	}

	void feed(const vector<double> &v, double alpha, double sigma) {
		//finding the closest point
		int ind_diff = -1;
		double diff = 1e9;
		for (int i = 0; i < cnt; ++i) {
			double d = dist(v, w[i], 1);
			if (diff > d) {
				diff = d;
				ind_diff = i;
			}
		}

		//the choice of neighbors
		//actually not yet

		//re-weighting
		for (int i = 0; i < cnt; ++i) {
			double k = alpha * exp(-dist(x[ind_diff], x[i], 0) / 2 / sigma / sigma);
			for (int j = 0; j < dim; ++j) {
				w[i][j] += k * (v[j] - w[i][j]);
			}
		}
	}

};




int main() {
	/*int tt = 1492556333;
	cout << tt << endl;*/
	srand(time(NULL));
	//srand(1492556333);


	string input_file_name = "229";
	freopen((R"(C:\Users\DNS\Desktop\)" + input_file_name + ".txt").c_str(), "r", stdin);

	int dim, w, h;
	cin >> dim >> h >> w;

	int cnt_epoh; cin >> cnt_epoh;

	double kalpha; cin >> kalpha;
	double ksigma; cin >> ksigma;


	int n;
	cin >> n;
	vector<vector<double> > a(n, vector<double>(dim, 0));
	vector<pair<double, double> > b(dim, make_pair(1000000, -1000000));
	vector<int> c(n);
	map<int, vector<double> > center;
	map<int, int> cnt_class;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < dim; ++j) {
			cin >> a[i][j];
			b[j].first = min(b[j].first, a[i][j]);
			b[j].second = max(b[j].second, a[i][j]);
		}
		cin >> c[i];
		auto it = center.find(c[i]);
		if (it == center.end()) {
			center[c[i]] = a[i];
			cnt_class[c[i]] = 1;
		}
		else {
			for (int j = 0; j < dim; ++j)
				it->second[j] += a[i][j];
			cnt_class[c[i]]++;
		}
	}

	for (auto it = center.begin(); it != center.end(); ++it) {
		for (int j = 0; j < dim; ++j) {
			it->second[j] /= cnt_class[it->first];
		}
	}
	h = sqrt(n);
	w = sqrt(n);
	vector<pair<int, int> > perm(n);
	for (int i = 0; i < n; ++i) {
		perm[i] = make_pair(rand(), i);
	}
	sort(perm.begin(), perm.end());


	kohonen_map kmp(dim, w, h, b);

	for (int i = 0; i < cnt_epoh; ++i) {
		if (i % 1000 == 0)
			cout << i << endl;
		int ind = i % n;
		double alpha = kalpha * (1 - (double)i / cnt_epoh) + 0.001;
		double sigma = ksigma * (1 - (double)i / cnt_epoh) + 2;
		kmp.feed(a[perm[ind].second], alpha, ksigma);
	}

#ifdef kmp_out
	ofstream fout;
	fout.open((R"(E:\Neyron\Koxonen_map\)" + input_file_name + "_out.txt").c_str());
	fout.precision(3);
	fout << fixed;

	for (int i = 0; i < kmp.w.size(); ++i) {
		for (int j = 0; j < dim; ++j) {
			fout << kmp.w[i][j] << ' ';
		}
		fout << 1;
		fout << endl;
	}
#endif

#ifdef kohon_map_out
	ofstream fkohon_out;
	fkohon_out.open((R"(E:\Neyron\Koxonen_map\)" + input_file_name + "_kohon_out.txt").c_str());
	fkohon_out.precision(3);
	fkohon_out << fixed;
	for (int i = 0; i < kmp.w.size(); ++i) {
		int ind_diff = -1, color = 0;
		double diff = 1e12;
		for (auto it = center.begin(); it != center.end(); ++it, ++color) {
			double d = dist(kmp.w[i], it->second, 0);
			if (diff > d) {
				diff = d;
				ind_diff = color;
			}
		}
		for (int j = 0; j < 2; ++j) {
			fkohon_out << kmp.x[i][j] << ' ';
		}
		fkohon_out << ind_diff << ' ' << ind_diff << endl;
	}

#endif

	return 0;
}
