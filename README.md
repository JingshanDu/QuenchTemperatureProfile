# QuenchTemperatureProfile

A script that estimates the heat transfer / cooling rate of a solid wafer/plate in a continuous binary mixed gas flow.

In this example, cooling of a silicon or glassy carbon plate in a hydrogen:argon mixture flow in a tube is calculated.

The resulting plot shows two curves: (1) plate temperature (assuming this is uniform across the plate), and (2) surface flow temperature (defined as the average of the plate and gas temperature).

## Notes

1. Heat transfer rates through different mechanisms are added up
2. Forced convection has a laminar flow
3. Data sources and calculation methods are referenced in each function
4. Solid conduction has not been implemented in this script

## Directions

Customize the gas flow and physical properties of the plates/gas to satisfy your system conditions. Enable or disable different heat transfer mechanisms in the Parameter Setup section. Then this code should run and plot the temperature profile.
