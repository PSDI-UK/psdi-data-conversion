# Conversion Chaining

## Motivation

The goal of this project is to provide a one-stop for conversion between the hundreds of file formats used in chemistry to describe molecular structure. This is done in part by using existing converters as plugins and allowing the user to pick a converter which can perform a desired conversion.

However, there are many pairs of formats where no direct conversion is possible, and we aim to handle this through chained conversions using multiple converters. This presents the problem: How do we determine possible chained conversions? And having determined these possible conversions, how do we determine the best one to use?

To answer this question, it helps to reframe the problem in terms of the mathematical structure known as the "graph", in particular the [directed graph](https://en.wikipedia.org/wiki/Directed_graph) (since there are some cases where conversion is only allowed in one direction), with formats represented as vertices and conversions as edges. Existing packages such as `igraph` allow us to take advantage of pre-existing code to solve problems such as this rather than needing to write our own from scratch.

The full graph of formats and allowed conversions looks like the following, with black dots representing different formats and different edge colors representing different converters performing each conversion:

![Graph of all conversions](img/all_conversions.png)

We can see from this that most formats are only supported by one converter, with only a handful being supported by more than one. This underlines the necessity of supporting chained conversions, as a random pair of two formats will not be supported by a direct conversion more often than not.

## Deciding an optimal path

Rather than tackling this full graph, let's look at a much smaller, more manageable subset it, with just four supported formats:

![Graph of four formats and their conversions, showing two bridging formats and two formats only supported by one converter](img/simple_graph.svg)

This graph shows four formats: MOLDY, PDB, CIF, and InCHi. PDB and CIF are supported by both the Atomsk (red) and Open Babel (blue) converters, while MOLDY is supported only by Atomsk and InCHi is supported only by Open Babel.

A user who wished to convert from MOLDY to InCHi would thus face a problem here, that no converter can perform this conversion directly. Our goal then is to determine a possible path the user could take. Looking at this graph, it's easy to intuit that there are two possible reasonable paths:

1. Convert from MOLDY to PDB with Atomsk, then from PDB to InCHi with Open Babel
2. Convert from MOLDY to CIF with Atomsk, then from CIF to InCHi with Open Babel

But with a much larger graph, it won't always be so easy to find a path. And in fact, even in this case, these are just the most obvious paths to a human. An exhaustive search would see many more paths, such as:

3. Convert from MOLDY to PDB with Atomsk, then from PDB to CIF with Atomsk, then from CIF to InCHi with Open Babel
4. Convert from MOLDY to CIF with Atomsk, then from CIF to PDB with Open Babel, then from PDB to InCHi with Open Babel

If we don't put in a constraint to avoid retracing one's steps, the number of paths is in fact infinite, but even without that constraint we see that a computer will pick up on many paths which are obviously suboptimal.

This is a common problem faced when working with graphs, with the most common real-world example of this coming up being pathfinding, i.e. finding directions from one location to another. Let's use this as an analogy to show how we would determine an optimal path.

We'll consider here four locations (analogous to the formats): New York, London, Edinburgh, and Oxford. There are flight connections between all of New York, London, and Edinburgh, and it's possible to drive between all of London, Edinburgh, and Oxford.

So our analogous pathfinding problem to the format conversion problem above would someone who wishes to travel from New York to Oxford. This can't be done solely by flying or by driving, but can be done with a hybrid of the two. However, which path is best depends on what the goal is. Let's say the goal is purely to minimise travel time. The travel times between each location are:

- New York and London (flight): 7 hours
- New York and Edinburgh (flight): 6 hours
- London and Edinburgh (flight): 1 hour
- London and Edinburgh (driving): 8 hours
- London and Oxford (driving): 2 hours
- Edinburgh and Oxford (driving): 7 hours

An exhaustive search of all paths that don't retrace their steps would find here that the fastest route is New York to London by flight (7 hours), then London to Oxford by driving (2 hours).

But an exhaustive search isn't reasonable with larger maps, as the number of possible routes scales exponentially with the number of vertices (`O(e^V)` time). Luckily, this is a solved problem in graph theory, and modern implementations of [Dijkstra's_algorithm](https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm) can find the optimal route in `O(E + V log(V))` time (where E is the number of edges and V is the number of vertices), making the problem easily tractable for a graph of the size we're working with.
