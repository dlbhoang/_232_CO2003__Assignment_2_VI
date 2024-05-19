#include "kDTree.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <queue>

kDTree::kDTree(int k) {
    this->k = k;
    root = nullptr;
}

kDTree::~kDTree() {
    destroyTree(root);
}

const kDTree &kDTree::operator=(const kDTree &other) {
    if (this != &other) {
        destroyTree(root);
        root = copyTree(other.root);
        k = other.k;
    }
    return *this;
}

kDTree::kDTree(const kDTree &other) {
    root = copyTree(other.root);
    k = other.k;
}

void kDTree::insert(const vector<int> &point) {
    root = insertRec(root, point, 0);
}

kDTreeNode* kDTree::insertRec(kDTreeNode *node, const vector<int> &point, int depth) {
    if (node == nullptr) {
        return new kDTreeNode(point);
    }

    int currentDimension = depth % k;

    if (point[currentDimension] < node->data[currentDimension]) {
        node->left = insertRec(node->left, point, depth + 1);
    } else {
        node->right = insertRec(node->right, point, depth + 1);
    }

    return node;
}

void kDTree::destroyTree(kDTreeNode *node) {
    if (node != nullptr) {
        destroyTree(node->left);
        destroyTree(node->right);
        delete node;
    }
}

kDTreeNode* kDTree::copyTree(kDTreeNode *node) {
    if (node == nullptr) {
        return nullptr;
    }

    auto newNode = new kDTreeNode(node->data);
    newNode->left = copyTree(node->left);
    newNode->right = copyTree(node->right);
    return newNode;
}

void kDTree::preorderTraversal() const {
    preorderTraversalHelper(root);
}

void kDTree::buildTree(const std::vector<std::vector<int>>& points) {
    std::vector<std::vector<int>> pointsCopy = points;
    buildTreeRec(root, pointsCopy, 0);
}

void kDTree::buildTreeRec(kDTreeNode*& node, std::vector<std::vector<int>>& points, int depth) {
    if (points.empty()) {
        node = nullptr;
        return;
    }

    int currentDimension = depth % k;

    std::sort(points.begin(), points.end(), [currentDimension](const std::vector<int>& a, const std::vector<int>& b) {
        return a[currentDimension] < b[currentDimension];
    });

    int medianIndex = points.size() / 2;
    node = new kDTreeNode(points[medianIndex]);

    std::vector<std::vector<int>> leftPoints(points.begin(), points.begin() + medianIndex);
    std::vector<std::vector<int>> rightPoints(points.begin() + medianIndex + 1, points.end());

    buildTreeRec(node->left, leftPoints, depth + 1);
    buildTreeRec(node->right, rightPoints, depth + 1);
}

bool kDTree::search(const vector<int>& point) const {
    return searchRec(root, point, 0);
}

bool kDTree::searchRec(kDTreeNode *node, const vector<int> &point, int depth) const {
    if (node == nullptr) {
        return false;
    }

    if (node->data == point) {
        return true;
    }

    int currentDimension = depth % k;

    if (point[currentDimension] < node->data[currentDimension]) {
        return searchRec(node->left, point, depth + 1);
    } else {
        return searchRec(node->right, point, depth + 1);
    }
}

kDTreeNode* kDTree::findMinRec(kDTreeNode* node, int d, int depth) {
    if (node == nullptr) {
        return nullptr;
    }

    int currentDimension = depth % k;

    if (currentDimension == d) {
        if (node->left == nullptr) {
            return node;
        }
        return findMinRec(node->left, d, depth + 1);
    }

    kDTreeNode* leftMin = findMinRec(node->left, d, depth + 1);
    kDTreeNode* rightMin = findMinRec(node->right, d, depth + 1);
    kDTreeNode* minNode = node;

    if (leftMin != nullptr && leftMin->data[d] < minNode->data[d]) {
        minNode = leftMin;
    }
    if (rightMin != nullptr && rightMin->data[d] < minNode->data[d]) {
        minNode = rightMin;
    }

    return minNode;
}

kDTreeNode* kDTree::removeRec(kDTreeNode* node, const vector<int>& point, int depth) {
    if (node == nullptr) {
        return nullptr;
    }

    int currentDimension = depth % k;

    if (node->data == point) {
        if (node->right != nullptr) {
            kDTreeNode* min = findMinRec(node->right, currentDimension, depth + 1);
            node->data = min->data;
            node->right = removeRec(node->right, min->data, depth + 1);
        } else if (node->left != nullptr) {
            kDTreeNode* min = findMinRec(node->left, currentDimension, depth + 1);
            node->data = min->data;
            node->right = removeRec(node->left, min->data, depth + 1);
            node->left = nullptr;
        } else {
            delete node;
            return nullptr;
        }
    } else if (point[currentDimension] < node->data[currentDimension]) {
        node->left = removeRec(node->left, point, depth + 1);
    } else {
        node->right = removeRec(node->right, point, depth + 1);
    }

    return node;
}

void kDTree::remove(const vector<int>& point) {
    root = removeRec(root, point, 0);
}

void kDTree::preorderTraversalHelper(kDTreeNode* node) const {
    if (node != nullptr) {
        std::cout << "(";
        for (int i = 0; i < node->data.size(); ++i) {
            std::cout << node->data[i];
            if (i < node->data.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << ") ";

        preorderTraversalHelper(node->left);
        preorderTraversalHelper(node->right);
    }
}


void kDTree::nearestNeighbourRec(kDTreeNode *node, const vector<int> &target, int depth, kDTreeNode *&best, double &bestDist) {
    if (node == nullptr) {
        return;
    }

    double dist = distance(node->data, target);
    if (dist < bestDist) {
        bestDist = dist;
        best = node;
    }

    int currentDimension = depth % k;
    kDTreeNode *nextNode = (target[currentDimension] < node->data[currentDimension]) ? node->left : node->right;
    kDTreeNode *otherNode = (target[currentDimension] < node->data[currentDimension]) ? node->right : node->left;

    nearestNeighbourRec(nextNode, target, depth + 1, best, bestDist);

    if (fabs(target[currentDimension] - node->data[currentDimension]) < bestDist) {
        nearestNeighbourRec(otherNode, target, depth + 1, best, bestDist);
    }
}

void kDTree::nearestNeighbour(const vector<int> &target, kDTreeNode *&best) {
    best = nullptr;
    double bestDist = std::numeric_limits<double>::max();
    nearestNeighbourRec(root, target, 0, best, bestDist);
}

void kDTree::kNearestNeighbourRec(kDTreeNode *node, const vector<int> &target, int depth, int k, priority_queue<pair<double, kDTreeNode *>> &pq) {
    if (node == nullptr) {
        return;
    }

    double dist = distance(node->data, target);
    if (pq.size() < k) {
        pq.push({dist, node});
    } else if (dist < pq.top().first) {
        pq.pop();
        pq.push({dist, node});
    }

    int currentDimension = depth % this->k;
    kDTreeNode *nextNode = (target[currentDimension] < node->data[currentDimension]) ? node->left : node->right;
    kDTreeNode *otherNode = (target[currentDimension] < node->data[currentDimension]) ? node->right : node->left;

    kNearestNeighbourRec(nextNode, target, depth + 1, k, pq);

    if (fabs(target[currentDimension] - node->data[currentDimension]) < pq.top().first || pq.size() < k) {
        kNearestNeighbourRec(otherNode, target, depth + 1, k, pq);
    }
}

void kDTree::kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList) {
    priority_queue<pair<double, kDTreeNode *>> pq;
    kNearestNeighbourRec(root, target, 0, k, pq);

    while (!pq.empty()) {
        bestList.push_back(pq.top().second);
        pq.pop();
    }

    reverse(bestList.begin(), bestList.end());
}

double kDTree::distance(const vector<int> &point1, const vector<int> &point2) const {
    double dist = 0;
    for (int i = 0; i < point1.size(); ++i) {
        dist += (point1[i] - point2[i]) * (point1[i] - point2[i]);
    }
    return sqrt(dist);
}



kNN::kNN(int k) {
    this->k = k;
}

void kNN::fit(Dataset &X_train, Dataset &y_train) {
    // Implementation of fit
}

Dataset kNN::predict(Dataset &X_test) {
    Dataset y_pred;
    // Implementation of predict
    return y_pred;
}

double kNN::score(const Dataset &y_test, const Dataset &y_pred) {
    double accuracy = 0.0;
    // Implementation of score
    return accuracy;
}