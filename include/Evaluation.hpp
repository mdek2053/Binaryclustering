#pragma once

#include <cmath>
#include <cstdint>
#include <set>
#include <vector>

// Contigency table
class CT {
public:
    std::vector<std::vector<int>> table;
    std::vector<int> rowsums;
    std::vector<int> colsums;
    CT(const std::vector<std::vector<int>>& table, const std::vector<int>& rowsums, const std::vector<int>& colsums) {
        this->table = table;
        this->rowsums = rowsums;
        this->colsums = colsums;
    }
};

// Confusion Matrix
class CM {
public:
    int64_t tp;
    int64_t fp;
    int64_t fn;
    int64_t tn;
    CM(int64_t tp, int64_t fp, int64_t fn, int64_t tn) {
        this->tp = tp;
        this->fp = fp;
        this->fn = fn;
        this->tn = tn;
    }
};

int countClusters(int nCells, const std::vector<int>& clusters) {
    std::set<int> ids;
    for (int i = 0; i < nCells; i++) {
        ids.insert(clusters[i]);
    }
    return ids.size();
}

int countClusters(const std::vector<int>& clusters) {
    return countClusters(clusters.size(), clusters);
}

double entropy(int n, const std::vector<int>& clusters) {
    if (n == 0) {
        return 1.0;
    }

    int nClusters = countClusters(n, clusters);

    if (nClusters == 1) {
        return 0.0;
    }

    auto counts = std::vector<int>(nClusters);
    for (int i : clusters) {
        counts[i]++;
    }

    double entropy = 0;
    for (int count : counts) {
        entropy -= (((double)count / n) * (log(count) - log(n)));
    }

    return entropy;
}

CT contengincyTable(int n, const std::vector<int>& clustersTrue, const std::vector<int>& clustersPred) {
    int nClustersTrue = countClusters(n, clustersTrue);
    int nClustersPred = countClusters(n, clustersPred);

    auto table = std::vector<std::vector<int>>(nClustersPred);
    for (int i = 0; i < nClustersPred; i++) {
        table[i] = std::vector<int>(nClustersTrue);
        for (int j = 0; j < nClustersTrue; j++) {
            int nij = 0;
            for (int k = 0; k < n; k++) {
                if (i == clustersPred[k] && j == clustersTrue[k]) {
                    nij++;
                }
            }
            table[i][j] = nij;
        }
    }

    auto rs = std::vector<int>(nClustersPred);
    auto cs = std::vector<int>(nClustersTrue);
    for (int i = 0; i < nClustersPred; i++) {
        for (int j = 0; j < nClustersTrue; j++) {
            rs[i] = rs[i] + table[i][j];
            cs[j] = cs[j] + table[i][j];
        }
    }

    return CT(table, rs, cs);
}

CM confusionMatrix(int n, const CT& table) {
    int64_t sumOfSquares = 0;
    for (const auto& v : table.table) {
        for (int i : v) {
            sumOfSquares += ((int64_t)i * i);
        }
    }

    int64_t tp = sumOfSquares - n;
    int64_t fp = -sumOfSquares;
    int64_t fn = -sumOfSquares;
    for (int i = 0; i < (int)table.rowsums.size(); i++) {
        for (int j = 0; j < (int)table.colsums.size(); j++) {
            fp += ((int64_t)table.table[i][j] * table.rowsums[i]);
            fn += ((int64_t)table.table[i][j] * table.colsums[j]);
        }
    }
    int64_t tn = (int64_t)n * n - fn - fp - sumOfSquares;
    return CM(tp, fp, fn, tn);
}

double calculateARI(CM cm) {
    if (cm.fp == 0 && cm.fn == 0) {
        return 1.0;
    }

    double numerator = (double)(cm.tp * cm.tn - cm.fn * cm.fp);
    double denominator = (double)((cm.tp + cm.fn) * (cm.fn + cm.tn) + (cm.tp + cm.fp) * (cm.fp + cm.tn));
    return 2.0 * numerator / denominator;
}

double calculateARI(int n, const std::vector<int>& clustersTrue, const std::vector<int>& clustersPred) {
    CT table = contengincyTable(n, clustersTrue, clustersPred);
    CM cm = confusionMatrix(n, table);
    return calculateARI(cm);
}

double calculateJaccard(CM cm) {
    return (double)cm.tp / (double)(cm.tp + cm.fn + cm.fp);
}

double calculateJaccard(int n, const std::vector<int>& clustersTrue, const std::vector<int>& clustersPred) {
    CT table = contengincyTable(n, clustersTrue, clustersPred);
    CM cm = confusionMatrix(n, table);
    return calculateJaccard(cm);
}

double calculateFM(CM cm) {
    return std::sqrt(((double)cm.tp / (double)(cm.tp + cm.fp)) * ((double)cm.tp / (double)(cm.tp + cm.fn)));
}

double calculateFM(int n, const std::vector<int>& clustersTrue, const std::vector<int>& clustersPred) {
    CT table = contengincyTable(n, clustersTrue, clustersPred);
    CM cm = confusionMatrix(n, table);
    return calculateFM(cm);
}

double calculateMI(int n, CT table, double entropyTrue) {
    double entropy = 0;
    for (int i = 0; i < (int)table.table.size(); i++) {
        std::vector<int> row = table.table[i];
        double rs = table.rowsums[i];
        double re = 0;

        for (int j = 0; j < (int)row.size(); j++) {
            if (table.table[i][j] != 0) {
                re += (table.table[i][j] / rs) * (log(table.table[i][j]) - log(rs));
            }
        }

        entropy -= ((rs / n) * re);
    }
    return entropyTrue - entropy;
}

double calculateMI(int n, const std::vector<int>& clustersTrue, const std::vector<int>& clustersPred) {
    CT table = contengincyTable(n, clustersTrue, clustersPred);
    double entropyTrue = entropy(n, clustersTrue);
    return calculateMI(n, table, entropyTrue);
}

double calculateNMI(int n, const std::vector<int>& clustersTrue, const std::vector<int>& clustersPred) {
    if (n == 0) {
        return 1.0;
    }

    double entropyTrue = entropy(n, clustersTrue);
    double entropyPred = entropy(n, clustersPred);

    if (entropyTrue == 0.0 || entropyPred == 0.0) {
        return 1.0;
    }

    CT table = contengincyTable(n, clustersTrue, clustersPred);
    double mi = calculateMI(n, table, entropyTrue);
    double homogeneity = mi / entropyTrue;
    double completeness = mi / entropyPred;

    if (homogeneity + completeness == 0.0) {
        return 0.0;
    } else {
        return (2.0 * homogeneity * completeness) / (homogeneity + completeness);
    }
}

/*
double calculateEMI(int n, CT table) {
        throw runtime_error("Not Implemented!");
}

double calculateAMI(int n, const vector<int>& clustersTrue, const vector<int>& clustersPred) {
        size_t nClustersTrue = countClusters(n, clustersTrue);
        size_t nClustersPred = countClusters(n, clustersPred);

        if ((nClustersPred == 1 && nClustersTrue == 1) || (nClustersPred == 0 && nClustersTrue == 0)) {
                return 1.0;
        }

        CT table = contengincyTable(n, clustersTrue, clustersPred);
        double hTrue = entropy(n, clustersTrue);
        double hPred = entropy(n, clustersPred);
        double mi = calculateMI(n, table, hTrue);
        double emi = calculateEMI(n, table);
        double normalizer = 0.5 * (hTrue + hPred);//mean
        double denominator = normalizer - emi;
        if (denominator < 0.0) {
                denominator = min(denominator, -numeric_limits<double>::epsilon());
        } else {
                denominator = max(denominator, numeric_limits<double>::epsilon());
        }
        return (mi - emi) / denominator;
}
*/

/*
double calculatePSI(int n, CT table) {
        int kr = table.rowsums.size();
        int kc = table.colsums.size();

        if (kr > kc) {
                throw runtime_error("CT needs to be transposed.");
        }

        vector<vector<double>> S = vector<vector<double>>(kr);
        for (int i = 0; i < kr; i++) {
                vector<int> row = table.table[i];
                S[i] = vector<double>(kc);

                for (int j = 0; j < kc; j++) {
                        S[i][j] = (double) row[j] / (double) max(table.rowsums[i], table.colsums[j]);
                        cout << S[i][j] << endl;
                }
        }

        vector<int> output_col4row2(kr);
        ssize_t retval = linear_sum_assignment(S.data(), xc, yc, output_col4row2.data(), false); // minimise=false

        double s = 0.0;
        for (ssize_t i = 0; i < xc; ++i)
                s += S[yc * i + output_col4row2[i]];

        std::sort(sum_x.begin(), sum_x.end());
        std::sort(sum_y.begin(), sum_y.end());
        double es = 0.0;
        for (ssize_t i = 0; i < xc; ++i)
                es += sum_y[yc - i - 1] * sum_x[xc - i - 1] / (double)std::max(sum_x[xc - i - 1], sum_y[yc - i - 1]);
        es /= (double)n;

        return max(0.0, (s - es) / (max(kr, kc) - es));
}
*/

/*
double calculatePSI(int n, const vector<int>& clustersTrue, const vector<int>& clustersPred) {
        CT table = contengincyTable(n, clustersTrue, clustersPred);
        return calculatePSI(n, table);
}
*/

/*
double calculateSilhouette(int nCells, const vector<int>& clusters, const vector<CellData>& data, SimilarityFunction metric) {
        size_t count = countClusters(nCells, clusters);
        double s = -1.0;
        if (count <= 1) {
                return s;
        }

        set<int> calculated;
        for (int i = 0; i < nCells; i++) {
                int cluster = clusters[i];
                vector<int> counts = vector<int>(count);
                vector<double> distances = vector<double>(count);
                if (calculated.find(cluster) == calculated.end()) {
                        for (int j = 0; j < nCells; j++) {
                                if (i != j) {
                                        int cluster2 = clusters[j];
                                        counts[cluster2]++;
                                        distances[cluster2] += metric.distance(data[i], data[j]);
                                }
                        }

                        double si = 0.0;
                        if (counts[cluster] > 0) {
                                double ai = distances[cluster] / counts[cluster];
                                double bi = numeric_limits<double>::max();
                                for (int j = 0; j < count; j++) {
                                        if (j != cluster) {
                                                bi = min(bi, distances[j] / counts[j]);
                                        }
                                }
                                si = (bi - ai) / max(ai, bi);
                        }

                        calculated.insert(cluster);
                        s = max(s, si);
                }
        }
        return s;
}
*/