# Conversion Chaining

## Motivation

The goal of this project is to provide a one-stop for conversion between the hundreds of file formats used in chemistry to describe molecular structure. This is done in part by using existing converters as plugins and allowing the user to pick a converter which can perform a desired conversion.

However, there are many pairs of formats where no direct conversion is possible, and we aim to handle this through chained conversions using multiple converters. This presents the problem: How do we determine possible chained conversions? And having determined these possible conversions, how do we determine the best one to use?

To answer this question, it helps to reframe the problem in terms of the mathematical structure known as the "graph", in particular the [directed graph](https://en.wikipedia.org/wiki/Directed_graph) (since there are some cases where conversion is only allowed in one direction - an example being the [InChIKey format](https://en.wikipedia.org/wiki/International_Chemical_Identifier#InChIKey), which is a hashed representation and thus can only be converted to, but not from), with formats represented as vertices and conversions as edges. Existing packages such as `igraph` allow us to take advantage of pre-existing code to solve problems such as this rather than needing to write our own from scratch.

The full graph of formats and allowed conversions looks like the following, with black dots representing different formats and different edge colors representing different converters performing each conversion:

![Graph of all conversions](img/all_conversions.png)

We can see from this that most formats are only supported by one converter, with only a handful being supported by more than one. This underlines the necessity of supporting chained conversions, as a random pair of two formats will not be supported by a direct conversion more often than not.

## General pathfinding

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

## Determining pathfinding weights

In the above example, we chose time as the relevant property we wanted to minimise in pathfinding, but of course we could have chosen something else, such as price, and minimised that instead. The chosen quantity is referred to as the "weight" of a path in the general case. Dijkstra's algorithm functions with any set of unbounded non-negative weights, even if they don't follow intuitive physical properties such as the [triangle inequality](https://en.wikipedia.org/wiki/Triangle_inequality), i.e. it's allowable for A->B to have a greater weight than the sum of A->C and C->B.

But what if we want multiple factors to be relevant, or what if there are other constraints we need to abide by, or wish to abide by if at all possible?

This is in fact the case we're dealing with when it comes to format conversions. Without getting into the details yet of what these mean, the goals we have are, in descending order of priority:

1. Don't follow a path which loses and then re-extrapolates any type of data shared by both the input and output format

2. Minimise loss of numerical accuracy (e.g. floating-point values being truncated in the conversion process)

3. Minimise time (implicitly minimising the number of steps)

Dijkstra's algorithm requires only one weight per path, so we have to find some way to combine these aspects.

### Weighted combination

The most straightforward way to have only one weight per path is to calculate it a weighted combination of the relevant factors, e.g.:

```math
\displaystyle W_i = w_{x}x_i + w_{y}y_i + w_{z}z_i
```

where $W_i$ is the total weight assigned to edge $i$ and $w_x$, $w_y$, and $w_z$ are the relative weights of factors $x$, $y$, and $z$ respectively. So we could for instance assign the highest weight to the most important factor, etc.

But the situation here is actually a bit more complicated. We don't simply wish for some of these factors to be weighted more, we wish for them to be strictly more important - no increase in factor $x$ can be compensated for by any decrease in $y$ or $z$. To use the travelling example, we could say that the high-importance factor is whether or not the traveller's luggage is permitted on the route, and the low-importance factor is how jostled the luggage will be in trip. Obviously if the luggage can't make the trip at all, it doesn't matter how smooth a ride it will have.

Strictly speaking, this can't be accomplished with a weighted combination, as no matter how different the weights are, there could always be extreme cases where the weight is overcome. This might occur rarely enough that it gives an acceptably low rate of error though, so we'll pin this possibility while we investigate if other solutions are worthwhile.

### Tiered pathfinding

In the example here of losing one's luggage, this is a binary event - either the condition of keeping one's luggage through the trip is satisfied or it isn't. With this binary condition, there isn't likely to be a single best path, but rather many paths which fulfill this criterion equally.

Pathfinding algorithms are capable of handling cases like this where there are many equally-good paths, so what we can do with this is, we can implement two stages of pathfinding. The first looks only at the high-importance criteria and identifies only the paths which satisfy it. The second stage looks at the second criteria and only the paths determined valid in the first step, and minimises the weight through these paths using the second criteria.

The fact that the high-importance criteria isn't binary actually isn't necessary for this solution. Instead of searching for all pathways which satisfy the condition, we can instead search for all pathways which minimise it to the same greatest extent.

This has the advantage over the previous approach that it will guarantee that the strict relative importance of the criteria is respective, but it comes with the drawback of greater computational overhead, needing to run the pathfinding algorithm multiple times (or else running some analogous operation on the list of shortest pathways from the first step such as a sort). This will also get more complicated to program if there are more than two importance tiers.

### Custom weight type

It's possible to run the pathfinding in a single stage with strict tiering of criteria if we use a custom data type for the weights of paths. The only requirements that Dijkstra's algorithm places on the weights is that they be non-negative, addable, and comparable. It's possible to construct a data type which meets these criteria and also allows for strict importance tiering, and in fact such a type is already in use for version numbering.

Version numbers are period-separated integers such as 0.1, 1.245.0, 0.2.40, etc. A difference in a more-major (earlier) number always takes precedence in a comparison over any difference in a less-major number, i.e. (X+1).0 is greater than X.Y for any value of Y, no matter how large, e.g. 1.0.0 is greater than 0.999999.0.

A number system such as this could be used for pathfinding with tiered importance simply by setting up appropriate weights, e.g.:

- Weight for luggage being allowed at all on the trip: 1.0.0
- Weight for time in hours of the trip: 0.2.0
- Weight for price in $100 for the trip: 0.1.0
- Weight for how jostled the luggage gets in the trip: 0.0.1

A single pathfinding algorithm could then be run, which will prioritise trips where luggage is allowed. Among those where it is (or among all if it isn't allowed on any route), it will balance time and price. If there are multiple best paths which tie on this as well, it will then prioritise whichever jostles the luggage the least.

This solution keeps the programming of the pathfinding simple (the extra complexity going into the definition of the data type), but will slow it down as comparisons of a custom data type such as this will take longer than native types, as compilers, hardware etc. are optimised for native numerical types. This also has the issue that if a third-party library is used for the pathfinding, it isn't likely to support a custom data type for weights. For instance, the `igraph` library only supports integer and floating-point weights.

## Optimal approach for our task

### Nature of the problem

To determine which approach is best for our task, let's now get into the details of what we need to do.

A chemical file format can store various types of information, of which we currently keep track of four, listed here in descending order of importance:

1. **Composition**: What elements make up the chemical (there are in fact some formats which don't store this information, such as [InChIKey](https://en.wikipedia.org/wiki/International_Chemical_Identifier#InChIKey))
2. **Connections**: Which atoms in the chemical are bonded to which other atoms
3. **Coordinates**: The physical locations of the atoms relative to each other. Some formats provide a full 3D information, while others only provide 2D information. To keep track of which is supported, our database lists 2D Coordinates and 3D Coordinates separately

At present, our database stores whether or not each of these properties is supported for many, but not all formats, listing the status as unknown for the remainder. While of course ideally this information should be added to the database, in the meantime it's an issue we need to take into account.

When a format is converted to another which is capable of storing more information than the source format, this information will either be excluded or extrapolated, depending on the specific scenario.

Different formats also store numeric information at different precisions (i.e. the number of digits), so a conversion to a format with lesser precision will result in some data loss. Even between formats with the same precision in theory, data loss may occur if it's represented differently (e.g. a point in 2D space could be described by its $x$ and $y$ coordinates or $r$ and $\theta$, and a conversion between the two will lose a small amount of information due to rounding at the final step), so it's safest to assume a small loss of information with every conversion.

All else being equal, the time to perform a conversion is also relevant.

For all considerations here aside from time, we also have to keep in mind whether the source and target format support this information or precision. For instance, consider three cases:

1. Format A (9-digit precision) -> Format B (9-digit precision) -> Format C (4-digit precision) -> Format D (9-digit precision)

2. Format E (4-digit precision) -> Format B (9-digit precision) -> Format C (4-digit precision) -> Format D (9-digit precision)

3. Format A (9-digit precision) -> Format B (9-digit precision) -> Format C (4-digit precision) -> Format F (4-digit precision)

Consider here the conversion from Format B to Format C, which goes from 9-digit precision to 4-digit precision. What should the weight be for the loss of data precision? In case 1, where the source format has 9-digit precision, this is the first step to lose that precision, and so should be weighted high to reflect this.

However, in case 2, the source format only has 4-digit precision, so although it's converted into a format with greater precision, it never actually has this precision to lose in the conversion from B to C, so it shouldn't be given a weight penalty for this. Case 3 shows the opposite issue: Although the source format does lose precision in this step, the target format doesn't support this precision, so the loss is inevitable and will have to happen at some point.

This leads to an important conclusion about this problem: Edge weights are not independent, instead being dependent on the source and target formats. In particular, we can say that a weight for losing information should be applied only if both the source and target format share that information (or level of precision).

If we want to be particularly careful, we should also say that this weight should only be applied the first time this information is lost along a pathway. This is difficult to implement though as it would mean that an edge would have a different weight depending on the path taken to get to its source vertex, which is not standard and would require a completely different pathfinding algorithm. However, this would require a particularly long chain of conversions - which is unlikely to occur - and the impact if it did would not be significant - an already dispreferred path would be even more dispreferred - so the amount of work needed to implement this would not nearly be justified.

### Constraint summary

Putting this all together, we have the following list of properties to weight, in descending order of importance:

1. Composition
2. Connections
3. 2D Coordinates
4. 3D coordinates
5. Numerical precision
6. Conversion time

With the following notes:

- We want to apply this importance strictly - e.g. a faster conversion isn't worth losing numerical precision
- We should always apply a small weight to numerical precision for each conversion to be on the safe side
- Not all formats currently list what information they support
- Any loss of information should only be weighted if both the source and target format support this information, and similarly with numerical precision

### Our approach

Let's look again at the graph of all conversions:

![Graph of all conversions](img/all_conversions.png)

Note that the vast majority of formats are supported by only one converter, with only a small number of formats being supported by more than one and acting as bridges. Converters also tend to support any-to-any conversions of formats they support, with the only exceptions being some formats they can only convert from, and others they can only convert to.

This implies that in the vast majority of cases, conversion chains will be short, likely with two intermediate formats at most, and typically just one. And since there are only a relative few formats that can act as bridges, any filter on equally-short pathways (however "short" is defined, as long as there's at least some cost to each conversion) is going to result in a relatively small number of pathways.

We also face the key issue here that any weights we assign to conversions (aside from for the conversion time) will depend on what the original source format and final target format are. This means there is going to be computational overhead in updating the edge weights, which will be more costly the more edges we have to deal with.

Since any initial filter will reduce the number of potential pathways so much, the additional computational overhead cost of the tiered pathfinding approach will be minimal. Meanwhile, with so many different factors to consider, the weighted combination approach could become unwieldy. Either solution is probably workable, but consideration of the factors involved leads us to the conclusion that a hybrid approach is best.

Due to the computational cost of calculating weights for all conversions, we want to do this on the full set of conversions only once, and then use this to reduce the number of formats and edges we need to worry about to a manageable number. After that point, it's less costly to make more complicated weight calculations.

Looking at the list the properties to weight, the first four factors (Composition, Connections, 2D Coordinates, and 3D coordinates) are all binary, the fifth (precision) requires a more complicated comparison and calculation for each conversion, and the sixth (time) is independent. This means we'll benefit most by cutting down the graph before dealing with precision. We thus apply a weighted combination approach first to the first four factors, cut down the graph, then use a second weighted combination approach for the final two.

### Weighted combination 1: Types of information

Since the first four factors are all binary, it can in fact be quite efficient to calculate weights for all of them at once by using bitwise operations. That is, we can use an unsigned integer where different bits represent whether or not we weight for the loss of a given property. We assign the following bits to these properties:

| Property       | Bit |
| -------------- | --- |
| Composition    | 24  |
| Connections    | 18  |
| 2D Coordinates | 12  |
| 3D Coordinates | 6   |
| (always)       | 0   |

Note here that we also include an always-weight bit at bit 0. This is to impose a minimal cost for any conversion, so that the pathfinding algorithm won't end up including circuitious paths through zero-weight conversions in the list of shortest paths returned from this step.

This set of bits was chosen so that they could all fit into a 32-bit integer value and they are spaced out enough that if used as weights there is negligible chance that any path will accrue enough weight from a lower-importance property to overcome the lack of a weight from a higher-importance property (the 6-bit difference would require 64 occurences to be overcome, which is exceedingly unlikely to cause a problem even imaging the future addition of many different specialised converters).

The procedure for this step is as follows:

1. For each conversion in the database, determine a bit mask, where 1 is set in the bit for a property if it's supported in the source property and not supported in the target property (if unknown for the source, assume it to be supported, and if unknown for the target, assume it to not be supported). This can be precalculated for every conversion, since this part doesn't depend on the source and target formats. For instance, take the conversion from molreport to InChI:

- Molreport supports Composition, Connections, and 2D Coordinates
- InChI supports Composition and Connections
- 2D Coordinates is thus supported only by the source format, so we set the bit for it (bit 12) to be 1, as well as the always-on bit (bit 0), getting the bit mask 00000000 00000000 00010000 00000001 (decimal representation 4,097)

2. Determine a bit mask for the specific conversion requested, where 1 is set in the bit for a property if it's supported in both the original source and final target format (if unknown, consider it supported). For instance, for the conversion from MMCIF to MOLDY:

- MMCIF supports all properties
- MOLDY has unknown support for all properties
- We assume MOLDY has support for all properties to be on the safe side. Thus, both the original source and final target support all properties, so we set 1 for all of them, getting the bit mask 00000001 00000100 00010000 01000001 (decimal representation 17,043,521)

3. Calculate the weights for all conversions by performing a bitwise "and" operation on their bit mask and the bit mask for the specific conversion requeted. This is the weight to be used for this conversion in determining a path for the conversion requested. For this example:

- The molreport to InChI conversion has bit mask 00000000 00000000 00010000 00000001
- The MMCIF to MOLDY conversion requested has bit mask 00000001 00000100 00010000 01000001
- A bitwise "and" combination of these gives 00000000 00000000 00010000 00000001 (decimal representation 4,097)

4. Run the pathfinding algorithm using these weights, collecting a list of the shortest paths

### Weighted combination 2: Precision and time

The first weighted combination will cut the number of possible shortest paths down to a much smaller number. In the current setup, there are likely to be fewer than 10 pathways at this point, so a full pathfinding approach might not even be necessary to determine the shortest based on loss of precision and time, with a simple search through the list potentially being faster given the lower overhead. However, this may not continue to be the case as the project expands and more converters are integrated with it, so a full pathfinding approach is probably still best here.

The first step here is to prune down the graph to just the formats and conversions which are on one of the shortest paths. This is best done by constructing sets of each from the list of shortest paths, then reconstructing a smaller graph using only these. Note that this smaller graph won't include all possible conversions between formats within it, but only those which were used in one of the shortest paths determined before.

Next, we have to determine a weight for precision and time for each of the edges in this new graph, prioritising the former. We'll start with the latter though, since placing an upper bound on it will let us know how much of a relative weight we need to give to precision so that differences in runtime will never dominate.

#### Time weights

At present, we lack information on the typical time to perform different conversions, so the following is speculative. While it is in theory possible to test and record times for all conversions, this is very demanding and likely not worth the efforts. Since all the converters incorporated so far support conversion between any pair of supported formats (albeit with some formats only being supported as sources and some only as targets), a rough formula to estimate the time for a given conversion would be:

```math
t_{\rm tot}(x) = t_{\rm o} + t_{\rm s}(x) + t_{\rm t}(x)
```

where $t_{\rm o}$ is the overhead time used by the converter for any conversion, $t_{\rm s}$ is the time needed to convert from the source format, $t_{\rm t}$ is the time needed to convert to the target format, and $x$ represents the size of the file to be converted.

The overhead time will be the same for all conversions, and we already intend for the precision weight to include a minimum weight for all conversions, so it serves no purpose including it here. The remaining terms will likely share the same time-dependence, so we can reasonably divide it out, getting a simplified time weight equation:

```math
w_{\rm tot} = w_{\rm s} + w_{\rm t}
```

where the time weight for a given conversion is simply the sum of time weights for the source and target formats (note that this will still differ for each converter). Reasonable values for these could be determined through the following procedure:

1. Choose a "standard" file. Convert it to every format supported by each converter for direct conversion many times, timing the conversion process and averaging for each target format. The lowest such time will be our zero-point; subtract it from all times. The remaining time in milliseconds can then be used as the target time weight for each format. Also note the average of these weights, which can be used as a default value.
2. For each of these converted files, convert it back to the original format of the standard file, again timing, averaging, and subtracting off the zero-point. This will give the source time weight for each format.
3. Repeat this procedure for each converter, getting separate source and target weights for each format with each converter.

Note that this process won't produce weights for the format of the original standard file. This will have to be determined in a separate step, using comparisons of conversions to/from it with other formats and their conversions to each other.

Milliseconds will likely be the most useful unit for time weights, based on preliminary testing.

Note that as time is the lowest-priority factor in pathfinding, gathering accurate time weights is a low-priority task. In the meantime, default weights can be used.

#### Precision weights

From the previous section, we noted that time weights will be measured in milliseconds. From testing, a reasonable upper bound for the time of a conversion is 10 seconds, but the weights will be measured on smaller files and thus likely to be much less than this. Since we only have two factors to worry about for this combination, we can be sure we're well in the clear by using a base precision weight of 65,536 (equivalent to a total conversion time of ~65.5 seconds on small test files, which should never be the case). For brevity, we will omit this factor from the discussion below.

When a number is expressed to a given number of digits, e.g. "3.510", the possible true number represented will be one that could round to it, here 3.5095 through 3.5105. With no other information, all values in this range are equally likely, so this can be described as a uniform distribution, which has a standard deviation proportional to the range spanned. In cases where this number is converted to another number with the same precision, this standard deviation can be used as a rough estimate of the loss of precision. It's possible that numbers will simply be copied over without any loss of precision, but it's best to be safest and impose a minimal cost for each conversion.

If this value is stored to fewer decimal places, e.g. for "3.51", the range becomes 3.505 through 3.515, increasing the standard deviation by a factor of 10 per decimal place lost. So if we wished to represent this directly in the weight, we could impose weight of $10^N$, where $N$ is the number of decimal places lost, with a minimum value of 1.

This level of accurately representing the variance isn't necessary, however. We're free to scale the values as we wish, as long as there isn't a risk of lower-tier changes overwhelming higher-tier changes. If we were to use a similar method to the weights for types of information and use bits for each level of precision lost, we would need to reserve 3 bits (a factor of 8) for each digit to be comfortable.

If we account for the values taken up by the time weight, this gives only enough room to represent 5 digits of precision loss in a 32-bit integer, which isn't likely to be enough - formats can range from just a few digits up to 16, or possibly more. It's safest to use a 64-bit integer here, giving us triple the room to store precision weight information, which should be enough for reasonable use cases.

#### Putting them together

Now, let's put these together into a specific procedure. The weight for this tier will be the sum of a weight for time (capped at a value of 65,535) and a weight for precision (setting bits from 2^16 = 65,536 upwards based on which digits are lost). The time weight is simply the sum of the calculated weights for the source and target formats, with the cap imposed. The precision weight is more complicated, calculated through the following procedure:

1. Determine the lower precision of the original source and final target formats - call this number of digits $D_{\rm 0}$
2. For each conversion, compare the precision of the target format $D_i$ to $D_{\rm 0}$, to get the number of digits lost $L_i = D_{\rm 0} - D_i$, constrained to the values $0 \leq L_i \leq 12$ - the minimum of 0 imposes a minimum cost for each conversion, and the maximum of 12 prevents it from overflowing a 64-bit integer.
3. The precision weight is then determined by setting bit $16 + 3L_i$ to 1, i.e. setting the value to $2^{16 + 3L_i}$

To illustrate this, let's imagine the following scenario, for the pathway A -> B -> C -> D:

- Format A: 8-digit precision, source time weight 10, target time weight 11
- Format B: 9-digit precision, source time weight 20, target time weight 25
- Format C: 4-digit precision, source time weight 5, target time weight 6
- Format D: 9-digit precision, source time weight 30, target time weight 35

The time weights for each conversion are:

- A -> B: Source time weight of A (10) + Target time weight of B (25) = 35
- B -> C: Source time weight of B (20) + Target time weight of C (6) = 26
- C -> D: Source time weight of C (5) + Target time weight of D (35) = 40

For the precision weight, $D_{\rm 0}$ is set to the lower precision of format A (8) and D (9): 8. The precision weights for each conversion are then calculated to be:

- A -> B: Target format's precision is $D_{\rm B}=9$; $D_{\rm 0} - D_{\rm B} = -1$; bounded to $L_{\rm B} = 0$; giving precision weight $2^{16+3*0}=65,536$
- B -> C: Target format's precision is $D_{\rm C}=4$; $D_{\rm 0} - D_{\rm C} = 4$; already within bounds, so $L_{\rm C} = 4$; giving precision weight $2^{16+3*4}=268,435,456$
- C -> D: Target format's precision is $D_{\rm D}=9$; $D_{\rm 0} - D_{\rm D} = -1$; bounded to $L_{\rm D} = 0$; giving precision weight $2^{16+3*0}=65,536$

The combined precision and time weights for each conversion are then:

- A -> B: Time weight (35) + Precision weight (65,536) = 65,571
- B -> C: Time weight (26) + Precision weight (268,435,456) = 268,435,482
- C -> D: Time weight (40) + Precision weight (65,536) = 65,576

The combined precision and time weights can then be assigned to all conversions in the reduced graph, and the pathfinding algorithm run on it to find all equally-shortest paths. In the case that one shortest path is found, it's returned as the best path. If more than one is found, whichever is the first in the list can be returned.
