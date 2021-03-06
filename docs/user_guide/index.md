
---
knit: "bookdown::render_book"
title: "AZRtools documentation"
author: ["Devin Pastoor", "Imad Mali"]
description: "AZRtools documentation"
site: bookdown::bookdown_site
documentclass: book
---

# User Guide {-}

The user guide documents the functionality for the [AZRsim](https://www.astrazeneca.com/) R package which can be used to simulate ordinary differential equations in C and NONMEM.

In order to simulate the model in the **AZRsim** package we have to,

1. Create an appropriate text file representation of the model.
2. Compile the model with the `create_model()` function.
3. Simulate the compiled model with the `simulate()` function.
