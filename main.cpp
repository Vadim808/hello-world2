#include<iostream>
#include<vector>
#include<algorithm>
#include<map>
#include<iomanip>

struct edge {
    int u, v;
    double time;
    edge(int a, int b, double c) {
        u = a;
        v = b;
        time = c;
    }
};

struct next_city {
    int city;
    double time;
    next_city(int c, double t) {
        city = c;
        time = t;
    }
};

std::vector<int> d;
std::vector<int> way;
static constexpr int INF = std::numeric_limits<int>::max();

void dijkstra(std::vector<std::vector<next_city>>& g, std::map<std::pair<int, int>, int>& po,
              std::vector<double>& p, int n, int k) {
    d.resize(n);
    for (int i = 1; i < n; ++i) {
        d[i] = INF;
    }
    std::vector<bool> was;
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
            if ((v != 0 && po[{v, g[v][j].city}] <= 0) ||
                (v == 0 && po[{v, g[v][j].city}] < k)) {
                int to = g[v][j].city;
                double len;
                if (po[{v, g[v][j].city}] < 0) {
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

std::vector<int> wass;

void dfs(std::vector<std::vector<next_city>>& g, std::map<std::pair<int, int>, int>& poi, int v, int K) {
    wass[v] = 1;
    if (v == K) {
        return;
    }
    for (size_t i = 0; i < g[v].size(); ++i) {
        if (wass[K] == 1) {
            return;
        }
        if (poi[{v, g[v][i].city}] == 1) {
            if (wass[K] == 1) {
                return;
            }
            way[g[v][i].city] = v;
            poi[{v, g[v][i].city}] = 0;
            poi[{g[v][i].city, v}] = 0;
            dfs(g, poi, g[v][i].city, K);
            if (wass[K] == 1) {
                return;
            }
        }
    }
}

void add_value3(std::map<std::pair<int, int>, int>& m, int c1, int c2, int c3, int v) {
    m[{c1, c3}] = v;
    m[{c3, c2}] = v;
    m[{c2, c3}] = v;
    m[{c3, c1}] = v;
}

void add_value2(std::map<std::pair<int, int>, int>& m, int c1, int c2, int v) {
    m[{c1, c2}] = v;
    m[{c2, c1}] = v;
}

void g_add_value(std::vector<std::vector<next_city>>& g, int c1, int c2, int c3,
                 next_city& nc1, next_city& nc2, next_city& nc3) {
    g[c3].push_back(nc1);
    g[c3].push_back(nc2);
    g[c1].push_back(nc3);
    g[c2].push_back(nc3);
}

int main() {
    std::cout.precision(15);
    int N, M, k;
    std::cin >> N >> M >> k;
    int m = M;
    int n = N;
    std::vector<std::vector<next_city>> g;
    std::vector<double> p;
    std::vector<edge> edges;
    std::map<std::pair<int, int>, int> points, ed, ed1;
    g.resize(n + 1);
    int city1, city2, city3;
    double time;
    next_city nc(1, 0);
    points[{1, 0}] = 0;
    points[{0, 1}] = 0;
    edge e(0, 1, 0);
    edges.push_back(e);
    g[0].push_back(nc);
    std::vector<std::vector<int>> nn;
    nn.resize(n + 1);
    for (int i = 0; i < n + 1; ++i) {
        nn[i].resize(n + 1);
    }
    for (int i = 0; i < M; ++i) {
        std::cin >> city1 >> city2 >> time;
        if (nn[city1][city2] == 1) {
            n += 1;
            city3 = n;
            add_value3(points, city1, city2, city3, 0);
            add_value3(ed, city1, city2, city3, i);
            m += 1;
            add_value3(ed1, city1, city2, city3, static_cast<int>(edges.size()));
            g.resize(n + 1);
            next_city nc1(city2, time / 2), nc2(city1, time / 2), nc3(city3, time / 2);
            g_add_value(g, city1, city2, city3, nc1, nc2, nc3);
            edge e1(city1, city3, time / 2), e2(city3, city2, time / 2);
            edges.push_back(e1);
            edges.push_back(e2);
        } else {
            nn[city1][city2] = 1;
            nn[city2][city1] = 1;
            add_value2(points, city1, city2, 0);
            add_value2(ed, city1, city2, i);
            add_value2(ed1, city1, city2, static_cast<int>(edges.size()));
            next_city nc1(city2, time), nc2(city1, time);
            edge e(city1, city2, time);
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
            std::pair<int, int> my_pair1(way[place], place), my_pair2(place, way[place]);
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
    if (points[{0, 1}] != k) {
        std::cout << -1;
    } else {
        double t = 0;
        std::vector<std::vector<int>> ways;
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
                ways[i].push_back(ed[{place, way[place]}]);
                t = t + edges[ed1[{place, way[place]}]].time;
                place = way[place];
            }
        }
        std::cout << t / k  << "\n";
        for (int i = 0; i < k; ++i) {
            size_t len = ways[i].size();
            for (size_t j = 0; j < ways[i].size(); ++j) {
                if (j > 0 && ways[i][j - 1] == ways[i][j]) {
                    len -= 1;
                    ways[i][j] = -1;
                }
            }
            std::cout << len << ' ';
            for (size_t j = 0; j < ways[i].size(); ++j) {
                if (ways[i][ways[i].size() - j - 1] != -1) {
                    std::cout << ways[i][ways[i].size() - j - 1] + 1 << ' ';
                }
            }
            std::cout << '\n';
        }
    }
}
