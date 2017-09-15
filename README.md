# crosslinked_monolayer
An mBuild + Foyer recipe for generating parameterized models of crosslinked alkylsilane monolayers attached to amorphous silica.

## Description
This recipe takes a three-stage approach to crosslinked monolayer creation:
    1. Adding and attaching surface-bound chains
    2. Adding crosslinked chains
    3. Identifying and adding crosslink bonds

Stages 1 & 2 utilize a steric approach to monolayer creation, under the assumption that monolayers
in experiment would pack to the maximum density without steric overlaps. For stage 1, this is
achieved by choosing available surface sites at random and attempting to attach chains.
A chain is only added to the system if no overlaps with existing chains would occur. Once all
remaining surface sites would yield overlaps with existing chains, the recipe enters stage 2.
Here, random positions within the *xy*-plane (i.e. the surface plane) are chosen for chain
insertion attempts. Again, chains are only inserted if no overlaps would occur with existing chains.
Once 10000 (or a different, user-defined number) consecutive failed insertion attempts are reached,
the monolayer is said to be at its maximum chain density.

After all monolayer chains have been added to the system, the crosslink network is determined.
This represents a problem similar to a "pickup-only" variant of the Multi-Depot Vehicle Routing Problem (MDVRP).
Here, the rules of our problem are that all chains must be traced through the crosslink
network back to a surface attachment site. The number of crosslink "clusters" or "strings"
is not defined. To solve this, the recipe first adds a crosslink between each chain and its
nearest neighbor. Following this, all subgraphs of the crosslink network are examined,
and subsequent crosslinks are created with nearest neighbors until all subgraphs feature
at least one chain attached directly to the surface.

## TO-DO
* Add description for how to download and install the package.
* Add a small code snippet to show off an example
* Add images of constructed monolayers and crosslink networks
* Add plot showing the fraction of surface bound chains as a function of spacing
