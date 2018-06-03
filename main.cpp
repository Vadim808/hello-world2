#include<iostream>
#include<vector>
#include<algorithm>
#include<map>
#include<iomanip>

using namespace std;

struct edge {
    int u, v;
    double time;
    void put(int a, int b, double c) {
        u = a;
        v = b;
        time = c;
    }
};

struct next_city {
    int city;
    double time;
    void put(int c, double t) {
        city = c;
        time = t;
    }
};

vector<int> d;
vector<int> way;
int INF = 1000000000;

void dijkstra(vector<vector<next_city>>& g, map<pair<int, int>, int>& po,
              vector<double>& p, int n, int k) {
    d.resize(n);
    for (int i = 1; i < n; ++i) {
        d[i] = INF;
    }
    vector<bool> was;
    was.resize(n);
    for (int i = 0; i < n; ++i) {
        int v = -1;
        for (int j = 0; j < n; ++j)
            if (!was[j] && (v == -1 || d[j] < d[v]))
                v = j;
        if (d[v] == INF)
            break;
        was[v] = true;
        
        for (size_t j = 0; j < g[v].size(); ++j) {
            if ((v != 0 && po[make_pair(v, g[v][j].city)] <= 0) ||
                (v == 0 && po[make_pair(v, g[v][j].city)] < k)) {
                int to = g[v][j].city;
                double len;
                if (po[make_pair(v, g[v][j].city)] < 0) {
                    len = -g[v][j].time + p[v] - p[to];
                } else {
                    len = g[v][j].time + p[v] - p[to];
                }
                if (d[v] + len < d[to]) {
                    d[to] = d[v] + len;
                    way[to] = v;
                }
            }
        }
    }
}

vector<int> wass;

void dfs(vector<vector<next_city>>& g, map<pair<int, int>, int>& poi, int v, int K) {
    wass[v] = 1;
    if (v == K) {
        return;
    }
    for (size_t i = 0; i < g[v].size(); ++i) {
        if (wass[K] == 1) {
            return;
        }
        if (poi[make_pair(v, g[v][i].city)] == 1) {
            if (wass[K] == 1) {
                return;
            }
            way[g[v][i].city] = v;
            poi[make_pair(v, g[v][i].city)] = 0;
            poi[make_pair(g[v][i].city, v)] = 0;
            dfs(g, poi, g[v][i].city, K);
            if (wass[K] == 1) {
                return;
            }
        }
    }
}

int main() {
    cout.precision(15);
    int N, M, k;
    cin >> N >> M >> k;
    int m = M;
    int n = N;
    vector<vector<next_city>> g;
    vector<double> p;
    vector<edge> edges;
    map<pair<int, int>, int> points;
    g.resize(n + 1);
    int city1, city2;
    double time;
    next_city nc;
    nc.put(1, 0);
    points[make_pair(1, 0)] = 0;
    points[make_pair(0, 1)] = 0;
    edge e;
    e.put(0, 1, 0);
    edges.push_back(e);
    g[0].push_back(nc);
    vector<vector<int>> nn;
    nn.resize(n + 1);
    for (int i = 0; i < n + 1; ++i) {
        nn[i].resize(n + 1);
    }
    map<pair<int, int>, int> ed;
    map<pair<int, int>, int> ed1;
    int city3;
    for (int i = 0; i < M; ++i) {
        cin >> city1 >> city2 >> time;
        if (nn[city1][city2] == 1) {
            n += 1;
            city3 = n;
            points[make_pair(city1, city3)] = 0;
            points[make_pair(city3, city2)] = 0;
            points[make_pair(city2, city3)] = 0;
            points[make_pair(city3, city1)] = 0;
            ed[make_pair(city1, city3)] = i;
            ed[make_pair(city3, city2)] = i;
            ed[make_pair(city3, city1)] = i;
            ed[make_pair(city2, city3)] = i;
            m += 1;
            ed1[make_pair(city1, city3)] = edges.size();
            ed1[make_pair(city3, city2)] = edges.size() ;
            ed1[make_pair(city3, city1)] = edges.size();
            ed1[make_pair(city2, city3)] = edges.size();
            g.resize(n + 1);
            next_city nc1, nc2, nc3;
            nc1.put(city2, time / 2); nc2.put(city1, time / 2);
            nc3.put(city3, time / 2);
            g[city3].push_back(nc1);
            g[city3].push_back(nc2);
            g[city1].push_back(nc3);
            g[city2].push_back(nc3);
            edge e1, e2;
            e1.put(city1, city3, time / 2);
            e2.put(city3, city2, time / 2);
            edges.push_back(e1);
            edges.push_back(e2);
        } else {
            nn[city1][city2] = 1;
            nn[city2][city1] = 1;
            points[make_pair(city1, city2)] = 0;
            points[make_pair(city2, city1)] = 0;
            ed[make_pair(city1, city2)] = i;
            ed[make_pair(city2, city1)] = i;
            ed1[make_pair(city1, city2)] = edges.size();
            ed1[make_pair(city2, city1)] = edges.size();
            next_city nc1, nc2;
            nc1.put(city2, time); nc2.put(city1, time);
            edge e;
            e.put(city1, city2, time);
            edges.push_back(e);
            g[city1].push_back(nc1);
            g[city2].push_back(nc2);
        }
    }
    p.resize(n + 1);
    for (int i = 0; i < n + 1; ++i) {
        p[i] = INF;
    }
    p[0] = 0;
    p[1] = 0;
    way.resize(n + 1);
    way[N] = -1;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m + 1; ++j) {
            if (p[edges[j].u] > p[edges[j].v] + edges[j].time) {
                way[edges[j].u] = edges[j].v;
                p[edges[j].u] = p[edges[j].v] + edges[j].time;
            }
            if (p[edges[j].v] > p[edges[j].u] + edges[j].time) {
                way[edges[j].v] = edges[j].u;
                p[edges[j].v] = p[edges[j].u] + edges[j].time;
            }
        }
    }
    int place;
    while (way[N] != -1) {
        place = N;
        while (true) {
            pair<int, int> my_pair1(way[place], place), my_pair2(place, way[place]);
            points[my_pair1] += 1;
            points[my_pair2] -= 1;
            if (points[my_pair1] == 1 && points[my_pair2] == 1) {
                points[my_pair1] = 0;
                points[my_pair2] = 0;
            }
            place = way[place];
            if (way[place] == 0) {
                my_pair1.first = 0;
                my_pair1.second = 1;
                points[my_pair1] += 1;
                break;
            }
        }
        dijkstra(g, points, p, n + 1, k);
        for (int i = 0; i < n + 1; ++i) {
            p[i] = p[i] + d[i];
        }
        way[N] = -1;
        dijkstra(g, points, p, n + 1, k);
    }
    if (points[make_pair(0, 1)] != k) {
        cout << -1;
    } else {
        double t = 0;
        vector<vector<int>> ways;
        ways.resize(k);
        
        for (int i = 0; i < k; ++i) {
            wass.clear();
            way.clear();
            way.resize(n + 1);
            wass.resize(n + 1);
            wass[0] = 1;
            dfs(g, points, 1, N);
            int place = N;
            while (way[place] != 0) {
                ways[i].push_back(ed[make_pair(place, way[place])]);
                t = t + edges[ed1[make_pair(place, way[place])]].time;
                place = way[place];
            }
        }
        cout << t / k  << "\n";
        for (int i = 0; i < k; ++i) {
            size_t len = ways[i].size();
            for (size_t j = 0; j < ways[i].size(); ++j) {
                if (j > 0 && ways[i][j - 1] == ways[i][j]) {
                    len -= 1;
                    ways[i][j] = -1;
                }
            }
            cout << len << ' ';
            for (size_t j = 0; j < ways[i].size(); ++j) {
                if (ways[i][ways[i].size() - j - 1] != -1) {
                    cout << ways[i][ways[i].size() - j - 1] + 1 << ' ';
                }
            }
            cout << '\n';
        }
    }
}
