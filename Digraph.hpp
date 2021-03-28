// Digraph.hpp
//
// ICS 46 Spring 2020
// Project #5: Rock and Roll Stops the Traffic
//
// This header file declares a class template called Digraph, which is
// intended to implement a generic directed graph.  The implementation
// uses the adjacency lists technique, so each vertex stores a linked
// list of its outgoing edges.
//
// Along with the Digraph class template is a class DigraphException
// and a couple of utility structs that aren't generally useful outside
// of this header file.
//
// In general, directed graphs are all the same, except in the sense
// that they store different kinds of information about each vertex and
// about each edge; these two types are the type parameters to the
// Digraph class template.

#ifndef DIGRAPH_HPP
#define DIGRAPH_HPP

#include <exception>
#include <functional>
#include <list>
#include <map>
#include <utility>
#include <vector>
#include <iostream>
#include <queue>


// DigraphExceptions are thrown from some of the member functions in the
// Digraph class template, so that exception is declared here, so it
// will be available to any code that includes this header file.

class DigraphException : public std::runtime_error
{
public:
    DigraphException(const std::string& reason);
};


inline DigraphException::DigraphException(const std::string& reason)
    : std::runtime_error{reason}
{
}



// A DigraphEdge lists a "from vertex" (the number of the vertex from which
// the edge points), a "to vertex" (the number of the vertex to which the
// edge points), and an EdgeInfo object.  Because different kinds of Digraphs
// store different kinds of edge information, DigraphEdge is a struct template.

template <typename EdgeInfo>
struct DigraphEdge
{
    int fromVertex;
    int toVertex;
    EdgeInfo einfo;
};



// A DigraphVertex includes two things: a VertexInfo object and a list of
// its outgoing edges.  Because different kinds of Digraphs store different
// kinds of vertex and edge information, DigraphVertex is a struct template.

template <typename VertexInfo, typename EdgeInfo>
struct DigraphVertex
{
    VertexInfo vinfo;
    std::list<DigraphEdge<EdgeInfo>> edges;
};



// Digraph is a class template that represents a directed graph implemented
// using adjacency lists.  It takes two type parameters:
//
// * VertexInfo, which specifies the kind of object stored for each vertex
// * EdgeInfo, which specifies the kind of object stored for each edge
//
// You'll need to implement the member functions declared here; each has a
// comment detailing how it is intended to work.
//
// Each vertex in a Digraph is identified uniquely by a "vertex number".
// Vertex numbers are not necessarily sequential and they are not necessarily
// zero- or one-based.

template <typename VertexInfo, typename EdgeInfo>
class Digraph
{
public:
    // The default constructor initializes a new, empty Digraph so that
    // contains no vertices and no edges.
    Digraph();

    // The copy constructor initializes a new Digraph to be a deep copy
    // of another one (i.e., any change to the copy will not affect the
    // original).
    Digraph(const Digraph& d);

    // The move constructor initializes a new Digraph from an expiring one.
    Digraph(Digraph&& d) noexcept;

    // The destructor deallocates any memory associated with the Digraph.
    ~Digraph() noexcept;

    // The assignment operator assigns the contents of the given Digraph
    // into "this" Digraph, with "this" Digraph becoming a separate, deep
    // copy of the contents of the given one (i.e., any change made to
    // "this" Digraph afterward will not affect the other).
    Digraph& operator=(const Digraph& d);

    // The move assignment operator assigns the contents of an expiring
    // Digraph into "this" Digraph.
    Digraph& operator=(Digraph&& d) noexcept;

    // vertices() returns a std::vector containing the vertex numbers of
    // every vertex in this Digraph.
    std::vector<int> vertices() const;

    // edges() returns a std::vector of std::pairs, in which each pair
    // contains the "from" and "to" vertex numbers of an edge in this
    // Digraph.  All edges are included in the std::vector.
    std::vector<std::pair<int, int>> edges() const;

    // This overload of edges() returns a std::vector of std::pairs, in
    // which each pair contains the "from" and "to" vertex numbers of an
    // edge in this Digraph.  Only edges outgoing from the given vertex
    // number are included in the std::vector.  If the given vertex does
    // not exist, a DigraphException is thrown instead.
    std::vector<std::pair<int, int>> edges(int vertex) const;

    // vertexInfo() returns the VertexInfo object belonging to the vertex
    // with the given vertex number.  If that vertex does not exist, a
    // DigraphException is thrown instead.
    VertexInfo vertexInfo(int vertex) const;

    // edgeInfo() returns the EdgeInfo object belonging to the edge
    // with the given "from" and "to" vertex numbers.  If either of those
    // vertices does not exist *or* if the edge does not exist, a
    // DigraphException is thrown instead.
    EdgeInfo edgeInfo(int fromVertex, int toVertex) const;

    // addVertex() adds a vertex to the Digraph with the given vertex
    // number and VertexInfo object.  If there is already a vertex in
    // the graph with the given vertex number, a DigraphException is
    // thrown instead.
    void addVertex(int vertex, const VertexInfo& vinfo);

    // addEdge() adds an edge to the Digraph pointing from the given
    // "from" vertex number to the given "to" vertex number, and
    // associates with the given EdgeInfo object with it.  If one
    // of the vertices does not exist *or* if the same edge is already
    // present in the graph, a DigraphException is thrown instead.
    void addEdge(int fromVertex, int toVertex, const EdgeInfo& einfo);

    // removeVertex() removes the vertex (and all of its incoming
    // and outgoing edges) with the given vertex number from the
    // Digraph.  If the vertex does not exist already, a DigraphException
    // is thrown instead.
    void removeVertex(int vertex);

    // removeEdge() removes the edge pointing from the given "from"
    // vertex number to the given "to" vertex number from the Digraph.
    // If either of these vertices does not exist *or* if the edge
    // is not already present in the graph, a DigraphException is
    // thrown instead.
    void removeEdge(int fromVertex, int toVertex);

    // vertexCount() returns the number of vertices in the graph.
    int vertexCount() const noexcept;

    // edgeCount() returns the total number of edges in the graph,
    // counting edges outgoing from all vertices.
    int edgeCount() const noexcept;

    // This overload of edgeCount() returns the number of edges in
    // the graph that are outgoing from the given vertex number.
    // If the given vertex does not exist, a DigraphException is
    // thrown instead.
    int edgeCount(int vertex) const;

    // isStronglyConnected() returns true if the Digraph is strongly
    // connected (i.e., every vertex is reachable from every other),
    // false otherwise.
    bool isStronglyConnected() const;

    // findShortestPaths() takes a start vertex number and a function
    // that takes an EdgeInfo object and determines an edge weight.
    // It uses Dijkstra's Shortest Path Algorithm to determine the
    // shortest paths from the start vertex to every other vertex
    // in the graph.  The result is returned as a std::map<int, int>
    // where the keys are vertex numbers and the value associated
    // with each key k is the precedessor of that vertex chosen by
    // the algorithm.  For any vertex without a predecessor (e.g.,
    // a vertex that was never reached, or the start vertex itself),
    // the value is simply a copy of the key.
    std::map<int, int> findShortestPaths(
        int startVertex,
        std::function<double(const EdgeInfo&)> edgeWeightFunc) const;


private:
    // Add whatever member variables you think you need here.  One
    // possibility is a std::map where the keys are vertex numbers
    // and the values are DigraphVertex<VertexInfo, EdgeInfo> objects.
    std::map<int, DigraphVertex<VertexInfo, EdgeInfo>> digraphMap;
    struct Vertex
    {
        int number;
        bool known = false;
        Vertex* predecessor = nullptr;
        double shortestDistance = std::numeric_limits<double>::infinity();
    };
    struct compare
    {
        bool operator () (Vertex v1, Vertex v2) const
        {
            return v1.shortestDistance < v2.shortestDistance;
        }
    };
    // You can also feel free to add any additional member functions
    // you'd like (public or private), so long as you don't remove or
    // change the signatures of the ones that already exist.
    bool checkVertexConnected(std::vector<int> visitedVertices, int vertex) const;
};



// You'll need to implement the member functions below.  There's enough
// code in place to make them compile, but they'll all need to do the
// correct thing instead.

template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph()
{
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph(const Digraph& d)
{
    digraphMap = d.digraphMap;
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph(Digraph&& d) noexcept
{
    std::swap(digraphMap, d.digraphMap);
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::~Digraph() noexcept
{
    digraphMap.clear();
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>& Digraph<VertexInfo, EdgeInfo>::operator=(const Digraph& d)
{
    if (this != &d)
    {
        digraphMap.clear();
        digraphMap = d.digraphMap;
    }
    return *this;
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>& Digraph<VertexInfo, EdgeInfo>::operator=(Digraph&& d) noexcept
{
    if (this != &d)
    {
        digraphMap.clear();
        std::swap(digraphMap, d.digraphMap);
    }
    return *this;
}


template <typename VertexInfo, typename EdgeInfo>
std::vector<int> Digraph<VertexInfo, EdgeInfo>::vertices() const
{
    std::vector<int> allVertices;
    for (typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::const_iterator 
    iter = digraphMap.begin(); iter != digraphMap.end(); iter++) //iterates through keys
    {
        allVertices.push_back(iter->first); //adds all keys to a vector
    }
    return allVertices;
}


template <typename VertexInfo, typename EdgeInfo>
std::vector<std::pair<int, int>> Digraph<VertexInfo, EdgeInfo>::edges() const
{
    std::vector<std::pair<int, int>> allEdges;
    for (typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::const_iterator 
    iter = digraphMap.begin(); iter != digraphMap.end(); iter++)
    {
        for (typename std::list<DigraphEdge<EdgeInfo>>::const_iterator iter2 = 
        iter->second.edges.begin(); iter2 != iter->second.edges.end(); iter2++)
        {
            allEdges.push_back(std::make_pair(iter2->fromVertex, iter2->toVertex));
            //get the edges from EdgeInfo and add them to a vector
        }
    }
    return allEdges;
}


template <typename VertexInfo, typename EdgeInfo>
std::vector<std::pair<int, int>> Digraph<VertexInfo, EdgeInfo>::edges(int vertex) const
{
    std::vector<std::pair<int, int>> allEdges;
    if (digraphMap.count(vertex) == 0)  //if vertex does not exist, throw exception
    {
        throw DigraphException("vertex does not exist");
    }
    //iterate through Digraph egdes for the specified vertex
    for (typename std::list<DigraphEdge<EdgeInfo>>::const_iterator iter = 
    digraphMap.at(vertex).edges.begin(); iter != digraphMap.at(vertex).edges.end(); iter++)
    {
            allEdges.push_back(std::make_pair(iter->fromVertex, iter->toVertex));
    }
    return allEdges;
}


template <typename VertexInfo, typename EdgeInfo>
VertexInfo Digraph<VertexInfo, EdgeInfo>::vertexInfo(int vertex) const
{
    if (digraphMap.count(vertex) == 0)  //if vertex does not exist, throw exception
    {
        throw DigraphException("vertex does not exist");
    }
    return digraphMap.at(vertex).vinfo;
}


template <typename VertexInfo, typename EdgeInfo>
EdgeInfo Digraph<VertexInfo, EdgeInfo>::edgeInfo(int fromVertex, int toVertex) const
{

    std::pair<int,int> targetEdge = std::make_pair(fromVertex,toVertex);
    std::vector<std::pair<int,int>> vertexEdges = edges(fromVertex); //gets all edges of fromVertex
    if (std::find(vertexEdges.begin(), vertexEdges.end(), targetEdge) == vertexEdges.end())
    {
        //checks for the (fromVertex, toVertex pair)
        throw DigraphException("edge does not exist");
    }
    for (typename std::list<DigraphEdge<EdgeInfo>>::const_iterator iter = 
    digraphMap.at(fromVertex).edges.begin(); iter != digraphMap.at(fromVertex).edges.end(); iter++)
    {
        if(iter->fromVertex == fromVertex && iter->toVertex == toVertex)
        {
            return iter->einfo;
        }
    }
    return EdgeInfo{};
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::addVertex(int vertex, const VertexInfo& vinfo)
{
    if (digraphMap.count(vertex) != 0)  //if vertex already exists, throw exception
    {
        throw DigraphException("vertex already exists");
    }
    DigraphVertex<VertexInfo, EdgeInfo> temp; //create digraph object
    temp.vinfo = vinfo;
    digraphMap.insert(std::pair<int,DigraphVertex<VertexInfo,EdgeInfo>>(vertex, temp));
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::addEdge(int fromVertex, int toVertex, const EdgeInfo& einfo)
{
    if (digraphMap.count(fromVertex) == 0 || digraphMap.count(toVertex) == 0)
    {
        throw DigraphException("vertex does not exist");
    }
    else 
    {
        for (typename std::list<DigraphEdge<EdgeInfo>>::const_iterator iter = 
        digraphMap.at(fromVertex).edges.begin(); iter != digraphMap.at(fromVertex).edges.end(); iter++)
        {
            if (iter->fromVertex == fromVertex && iter->toVertex == toVertex)
            {
                throw DigraphException("edge already exists");
            }
        }
    }
    DigraphEdge<EdgeInfo> someEdge;
    someEdge.fromVertex = fromVertex;
    someEdge.toVertex = toVertex;
    someEdge.einfo = einfo;
    digraphMap.at(fromVertex).edges.push_back(someEdge);
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::removeVertex(int vertex)
{
    if (digraphMap.count(vertex) == 0)  //if vertex not in digraphMap, throw exception
    {
        throw DigraphException("vertex does not exist");
    }
    for (typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::iterator 
    iter = digraphMap.begin(); iter != digraphMap.end(); iter++)
    {
        for (typename std::list<DigraphEdge<EdgeInfo>>::iterator iter2 = 
        iter->second.edges.begin(); iter2 != iter->second.edges.end(); iter2++)
        {
            if (iter2->toVertex == vertex) //find matching vertex
            {
                iter->second.edges.erase(iter2); //delete the edges of that vertex
            }
        }
    }
    digraphMap.erase(vertex); //delete the vertex
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::removeEdge(int fromVertex, int toVertex)
{
    if (digraphMap.count(fromVertex) == 0 || digraphMap.count(toVertex) == 0)
    {
        throw DigraphException("vertex does not exist");  //if vertex not in DigraphMap, throw exception
    } 
    for (typename std::list<DigraphEdge<EdgeInfo>>::iterator iter = 
    digraphMap.at(fromVertex).edges.begin(); iter != digraphMap.at(fromVertex).edges.end(); iter++)
    {
        if (iter->fromVertex == fromVertex && iter->toVertex == toVertex)
        {
            digraphMap.at(fromVertex).edges.erase(iter);
        }
    }
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::vertexCount() const noexcept
{
    return digraphMap.size();
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::edgeCount() const noexcept
{
    int totalEdges = 0;
    //iterate through each vertex and add all edges to totalEdges
    for (typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::const_iterator 
    iter = digraphMap.begin(); iter != digraphMap.end(); iter++)
    {
        totalEdges += iter->second.edges.size();
    }
    return totalEdges;
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::edgeCount(int vertex) const
{
    if (digraphMap.count(vertex) == 0)  //if vertex not in digraphMap, throw exception
    {
        throw DigraphException("vertex does not exist");
    }
    return digraphMap.at(vertex).edges.size();
}


template <typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::isStronglyConnected() const
{
    std::vector<int> allVertices = vertices();
    for (int i = 0; i < allVertices.size(); i++) //pass in each vertex
    {
        std::vector<int> visitedVertices;
        if (checkVertexConnected(visitedVertices, allVertices[i]) == false)
        {
            return false; //breaks if one vertex isn't connected
        }
    }
    return true;
}

template <typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::checkVertexConnected(
    std::vector<int> visitedVertices, int vertex) const
{
    std::vector<int> allVertices = vertices();
    for (typename std::list<DigraphEdge<EdgeInfo>>::const_iterator iter = 
    digraphMap.at(vertex).edges.begin(); iter != digraphMap.at(vertex).edges.end(); iter++)
    {
        //iterate through all edges of a vertex
        if (std::find(visitedVertices.begin(), visitedVertices.end(), 
        iter->toVertex) == visitedVertices.end()) 
        //check if the toVertex is not already in the vector
        {
            visitedVertices.push_back(iter->toVertex);
        }
        else 
        {
            checkVertexConnected(visitedVertices, iter->toVertex);
            
        }
    }
    if (digraphMap.size() == visitedVertices.size())
    {
        return true;
    }
    return false;
}

template <typename VertexInfo, typename EdgeInfo>
std::map<int, int> Digraph<VertexInfo, EdgeInfo>::findShortestPaths(
    int startVertex,
    std::function<double(const EdgeInfo&)> edgeWeightFunc) const
{
    std::map<int, int> shortestPaths;
    std::map<int, Vertex> allVertices;
    for (typename std::map<int,DigraphVertex<VertexInfo,EdgeInfo>>::const_iterator 
    iter = digraphMap.begin(); iter!= digraphMap.end(); iter++)
    {
        Vertex someVertex;
        someVertex.number = iter->first;
        allVertices.insert(std::pair<int, Vertex>(iter->first, someVertex));
    }
    Vertex firstVertex;
    firstVertex.number = startVertex;
    firstVertex.shortestDistance = 0;
    shortestPaths[startVertex] = startVertex;
    std::priority_queue <Vertex, std::vector<Vertex>, compare> pq;
    pq.push(firstVertex);
    while (pq.empty() == false)
    {
        Vertex smallestPriorityVertex = pq.top();
        pq.pop(); //dequeue smallestPriorityVertex
        if (smallestPriorityVertex.known == false)
        {
            smallestPriorityVertex.known = true;
            for(typename std::list<DigraphEdge<EdgeInfo>>::const_iterator iter = 
            digraphMap.at(smallestPriorityVertex.number).edges.begin();
            iter != digraphMap.at(smallestPriorityVertex.number).edges.end(); iter++)
            {
                if (allVertices.at(iter->toVertex).shortestDistance >
                firstVertex.shortestDistance + edgeWeightFunc(iter->einfo))
                {
                    shortestPaths[iter->toVertex] = smallestPriorityVertex.number;
                    allVertices.at(iter->toVertex).shortestDistance =
                    firstVertex.shortestDistance + edgeWeightFunc(iter->einfo);
                    allVertices.at(iter->toVertex).predecessor = &smallestPriorityVertex;
                    pq.push(allVertices.at(iter->toVertex));
                }
            }
        }
    }
    return shortestPaths;
}



#endif // DIGRAPH_HPP

