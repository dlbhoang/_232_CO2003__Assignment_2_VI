#ifndef KDTREE_HPP
#define KDTREE_HPP

#include <vector>
#include <iostream>
#include <queue>
#include "Dataset.hpp"

using namespace std;

struct kDTreeNode {
    vector<int> data;
    kDTreeNode *left;
    kDTreeNode *right;
    
    kDTreeNode(vector<int> data, kDTreeNode *left = nullptr, kDTreeNode *right = nullptr) {
        this->data = data;
        this->left = left;
        this->right = right;
    }

    friend ostream &operator<<(ostream &os, const kDTreeNode &node) {
        os << "(";
        for (int i = 0; i < node.data.size(); i++) {
            os << node.data[i];
            if (i != node.data.size() - 1) {
                os << ", ";
            }
        }
        os << ")";
        return os;
    }
};

class kDTree {
private:
    int k;
    kDTreeNode *root;

    kDTreeNode* insertRec(kDTreeNode *node, const vector<int> &point, int depth);
    kDTreeNode* removeRec(kDTreeNode *node, const vector<int> &point, int depth);
    kDTreeNode* findMinRec(kDTreeNode *node, int d, int depth);
    void destroyTree(kDTreeNode *node);
    kDTreeNode* copyTree(kDTreeNode *node);
    void nearestNeighbourRec(kDTreeNode *node, const vector<int> &target, int depth, kDTreeNode *&best, double &bestDist);
    void inorderTraversalRec(kDTreeNode *node) const;
    void preorderTraversalRec(kDTreeNode *node) const;
    void postorderTraversalRec(kDTreeNode *node) const;
    int heightRec(kDTreeNode *node) const;
    int nodeCountRec(kDTreeNode *node) const;
    int leafCountRec(kDTreeNode *node) const;
    bool searchRec(kDTreeNode *node, const vector<int> &point, int depth) const;
    double distance(const vector<int> &point1, const vector<int> &point2) const;
    void kNearestNeighbourRec(kDTreeNode *node, const vector<int> &target, int depth, int k, priority_queue<pair<double, kDTreeNode *>> &pq);

public:
    kDTree(int k = 2);
    ~kDTree();

    const kDTree &operator=(const kDTree &other);
    kDTree(const kDTree &other);

    void inorderTraversal() const;
    void preorderTraversal() const;
    void postorderTraversal() const;
    int height() const;
    int nodeCount() const;
    int leafCount() const;

    void insert(const vector<int> &point);
    void remove(const vector<int> &point); // To be implemented
    bool search(const vector<int> &point) const;
    void buildTree(const vector<vector<int>> &pointList);
    void nearestNeighbour(const vector<int> &target, kDTreeNode *&best);
    void kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList);
    vector<int> findNearest(const vector<int> &point);
    void preorderTraversalHelper(kDTreeNode* node) const;
    void buildTreeRec(kDTreeNode*& node, std::vector<std::vector<int>>& points, int depth);
};

class kNN {
private:
    int k;

public:
    kNN(int k = 5);

    void fit(Dataset &X_train, Dataset &y_train);
    Dataset predict(Dataset &X_test);
    double score(const Dataset &y_test, const Dataset &y_pred);

private:
    double distance(const vector<int> &point1, const vector<int> &point2) const;
};

#endif // KDTREE_HPP
