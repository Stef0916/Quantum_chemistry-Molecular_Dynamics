# Analysis of Relationship Between Variables

This document provides an overview of the relationship between the variables ΘΘ, αα, and JJ. It discusses the feature engineering process, the derived equation, and the nature of solutions we observed.

## 1. Feature Engineering

In order to capture the nuances of how Theta and Alpha relate to J, we employed a polynomial feature expansion. Here's a breakdown of the features:

- **X and Y**: Represents Theta and Alpha respectively.
- **X^2 and Y^2**: Squared terms of Theta and Alpha to capture potential non-linear relationships.
- **X*Y**: Interaction term to capture synergistic or antagonistic effects.
- **X^2*Y and X*Y^2**: Represents interaction effects where one variable is squared.
- **X^2*Y^2**: Combines squared effects of both variables.

## 2. Derived Equation

The polynomial regression resulted in the following equation:

J = 7498.5893663009 - 97.8261709408X - 74.8697514292Y + 0.1807431182X^2 - 0.0063908009X^2Y - 0.0003704683X^2Y^2 - 4.2194914164Y^2 + 0.0794858425XY^2 + 1.4614280743XY

## 3. Solution for the Quadratic Equation

The last step is to solve the equation to extract Y (= Alpha) given experimental values of X (= Theta) and J. These new values of Alpha were then used to compute new J values to asses the predictive power of the correlation.

## 4. Complete project

Check the project here: https://stefani-gamboa-portfolio.netlify.app/project2

## Contributing

This work was performed under the supervision of Dr. Mylis Orio within the framework of my Ph.D thesis in the Aix-Marseille University.
