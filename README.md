# Julia Life Cycle Model
We simulate a simple life cycle model in Julia. This is based on Matlab code from the course Dynamic Economics in Practice by Costa Dias and O'Dea. https://www.ifs.org.uk/publications/13801 We replicate models versions 1 and 5 from Costa Dias and O'Dea.

## Model v1
This program solves and simulates a finite period consumption and saving (cake-eating problem) problem. There is no income or uncertainty. The consumer starts with an endowment of assets and chooses how to allocate consumption over time. Solution is by value function iteration.

## Model v5
This version adds a realistic income stream with uncertainty. The household gradually accumulates assets to smooth consumption over income fluctuations and during retirement. 

# Useful Julia resources include:
* https://en.wikibooks.org/wiki/Introducing_Julia
* https://juliadocs.github.io/Julia-Cheat-Sheet/
* https://docs.julialang.org/

The julia code was written by Patrick Moran and David Sturrock
* https://sites.google.com/view/patrickmoran/
* https://sites.google.com/site/davidjsturrock/home
