#pragma once

#include <igraph/igraph.h>
#include <igraph/igraph_community.h>

#include <set>
#include <stdexcept>
#include <vector>

#include "Graph.hpp"

enum CommunityDetection {
    Louvain,
    Leiden,
};

std::string getName(CommunityDetection algorithm) {
    switch (algorithm) {
        case Louvain:
            return "Louvain";
        case Leiden:
            return "Leiden";
        default:
            throw std::runtime_error("Unknown community detection algorithm");
    }
}

CommunityDetection communityFromName(const std::string& name) {
    if (name == "Louvain") return Louvain;
    if (name == "Leiden") return Leiden;

    throw std::runtime_error("Unknown name for community detection algorithm");
}

std::vector<int> louvain(const Graph& g, float res) {
    /* Select edges */
    std::set<uint64_t> ids;
    std::vector<Edge> edgeList;
    for (int j = 0; j < g.size; j++) {
        for (Edge e : g.nodes[j].out) {
            uint64_t id = e.undirectedId();
            if (ids.find(id) == ids.end()) {
                edgeList.push_back(e);
                ids.insert(id);
            }
        }
    }
    ids.clear();

    /* Create edges */
    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, edgeList.size() * 2);
    igraph_vector_t weights;
    igraph_vector_init(&weights, edgeList.size());
    for (int i = 0; i < (int)edgeList.size(); i++) {
        VECTOR(edges)[i * 2] = edgeList[i].from;
        VECTOR(edges)[i * 2 + 1] = edgeList[i].to;
        VECTOR(weights)[i] = edgeList[i].weight;
    }
    edgeList.clear();

    /* Create graph */
    igraph_t graph;
    igraph_create(&graph, &edges, g.size, IGRAPH_UNDIRECTED);
    igraph_vector_int_destroy(&edges);

    /* Output */
    igraph_vector_int_t membership;
    igraph_vector_int_init(&membership, 0);

    /* Run Louvain */
    igraph_community_multilevel(&graph, &weights, res, &membership, NULL, NULL);

    /* Parse results */
    auto result = std::vector<int>(g.size);
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(&membership); i++) {
        result[i] = (int)VECTOR(membership)[i];
    }

    /* Destroy */
    igraph_destroy(&graph);
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&membership);

    /* Result */
    return result;
}

std::vector<int> leiden(const Graph& g, float res, float beta, int iterations) {
    /* Select edges */
    std::set<uint64_t> ids;
    std::vector<Edge> edgeList;
    for (int j = 0; j < g.size; j++) {
        for (Edge e : g.nodes[j].out) {
            uint64_t id = e.undirectedId();
            if (ids.find(id) == ids.end()) {
                edgeList.push_back(e);
                ids.insert(id);
            }
        }
    }
    ids.clear();

    /* Create edges */
    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, edgeList.size() * 2);
    igraph_vector_t weights;
    igraph_vector_init(&weights, edgeList.size());
    for (int i = 0; i < (int)edgeList.size(); i++) {
        VECTOR(edges)[i * 2] = edgeList[i].from;
        VECTOR(edges)[i * 2 + 1] = edgeList[i].to;
        VECTOR(weights)[i] = edgeList[i].weight;
    }
    edgeList.clear();

    /* Create graph */
    igraph_t graph;
    igraph_create(&graph, &edges, g.size, IGRAPH_UNDIRECTED);
    igraph_vector_int_destroy(&edges);

    /* Output */
    igraph_vector_int_t membership;
    igraph_vector_int_init(&membership, 0);

    /* Run Leiden */
    igraph_community_leiden(&graph, &weights, NULL, res, beta, false, iterations, &membership, NULL, NULL);

    /* Parse results */
    auto result = std::vector<int>(g.size);
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(&membership); i++) {
        result[i] = (int)VECTOR(membership)[i];
    }

    /* Destroy */
    igraph_destroy(&graph);
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&membership);

    /* Result */
    return result;
}

std::vector<int> clusterize(const Graph& g, CommunityDetection algorithm, float res, float beta = 0.01f, int iterations = -1) {
    switch (algorithm) {
        case Louvain:
            return louvain(g, res);
        case Leiden:
            return leiden(g, res, beta, iterations);
        default:
            throw std::runtime_error("Community detection algorithm unknown...");
    }
}
