
# AZRnlme 

## Running a model

A full workflow may look like the following:

```
# Load the create_model
model         <- create_model("model.txt")

# Define type of dosing that appears in the model
dosing        <- list(
  INPUT1 = c(type="BOLUS")
)

# Define path for NONMEM model export
projectPath   <- "MODEL_01"

# Define data information
data          <- list(
  relPathFromProject = "..",
  fileName           = "data.csv",
  headerIdent        = "IGNORE,ID,IGNORE,IGNORE,IGNORE,TIME,TIMEPOS,TAD,IGNORE,IGNORE,YTYPE,IGNORE,DV,IGNORE,CENS,MDV,EVID,AMT,ADM,IGNORE,RATE,IGNORE,IGNORE,COV,CAT"
)

# Define model specifications related to the parameter, statistical, covariate model
modelSpec     <- list(
  # Fixed effect
  POPvalues0     = c(ka=2, CL=40, Vc=2),
  POPestimate    = c(ka=1, CL=1,  Vc=1),

  # Random effects
  IIVvalues0     = c(ka=0.3, CL=0.4, Vc=0.5),
  IIVestimate    = c(ka=1, CL=1,  Vc=1),
  IIVdistribution = c(ka="L", CL="L", Vc="L")
)

# Create the NONMEM model
create_project(model,
                       dosing,
                       projectPath,
                       data,
                       modelSpec             = modelSpec)
# Run the NONMEM model
run_project(projectPath)
```

also of note, is that parameters much be declared as
being an `estimate` or a `regressor`

```
********** MODEL PARAMETERS
ka = 1   % (1/hours) <estimate>
CL = 3   % (L/hours) <estimate>
Vc = 60  % <estimate>
```

## dosing

Dosing can be specified of the type:

* BOLUS
* INFUSION
* ABSORPTION0

and is linked to an input

```
dosing        <- list(
  INPUT1 = c(type="BOLUS")
)
```

## project path

The project path sets where the nonmem model and data should be placed for execution.
If the path does not exist, it will be created. If it does exist, the files will be overwritten
when a model is created

It may be stored as a variable (recommended) or can be directly referenced.

```
projectPath   <- "MODEL_01"
```

## data

Define data information

There are a number of data requirements that must be met:

* IXGDF - Index column
* USUBJID - Subject ID
* ID - nonmem-style ID
* TIME - Time
* TIMEPOS - Time, with no times less than 0
* TAD - Time after dose
* DV - primary dependent variable
* MDV - missing dependent variable flag
* EVID - nonmem event ID
* CENS - censoring flag
* AMT - amount
* ADM - administration flag
* RATE - rate flag

Given a dataset, the location of these columns, and any others 
to be kept or ignored (replace name with IGNORE) must be provided 
in the `data` object as a comma separated list with no spaces called
`headerIndent`.

```
data          <- list(
  relPathFromProject = "..",
  fileName           = "data.csv",
  headerIdent        = "IGNORE,ID,IGNORE,IGNORE,IGNORE,TIME,TIMEPOS,TAD,IGNORE,IGNORE,YTYPE,IGNORE,DV,IGNORE,CENS,MDV,EVID,AMT,ADM,IGNORE,RATE,IGNORE,IGNORE,COV,CAT"
)
```

## model specification:

* POPestimate
* POPvalues0
* IIVdistribution
* IIVestimate
* IIVvalues0
* errorModel
* covarianceModel
* covariateModel
* covariateModelValues
* COVestimate
* COVcentering

with less common arguments also accepted:

* METHOD
* MAXEVAL
* SIGDIGITS
* PRINT
* M4
* SEED
* K1
* K2
* NRCHAINS
* IMPORTANCESAMPLING
* IMP_ITERATIONS
* ITS
* ITS_ITERATIONS

The common details are described below

### `POPvalues0` and `POPestimate`

POPvalues and estimate set the initial values and whether to estimate

```
modelSpec     <- list(
  # ... 
  POPvalues0     = c(ka=2, CL=40, Vc=2),
  POPestimate    = c(ka=1, CL=1,  Vc=1),
)
```


### `IIVdistribution`

because the underlying ODE system does not support heirarchical constructs,
the value of a single parameter can be overloaded to treat it as a 
heirarchical parameter by declaring one of 3 supported distribution,
`L`ognormal, `N`ormal, and `G`-logitnormal


### `IIVvalues0` and  `IIVestimate`

The IIVestimate specifies the initial estimate of the IIV parameters

```
modelSpec     <- list(
  # ... 
  IIVvalues0 = c(ka=1, CL =1.4, Vc = 10.8),
  IIVestimate    = c(ka=1, CL=1,  Vc=1)
)
```

### `errorModel`

Error models are linked to outputs and can be either

`abs`, `rel`, or `absrel` and must be specified as a list.

```
modelSpec <- list(
  # ...
  errorModel = list(
    OUTPUT1 = c("absrel", c(abs=1,rel=0.1))
  )
)
```

### `covarianceModel`

Given a need to set some components of the covariance model
to be a block, and others to be diagnoal, the block elements 
may be passed as a vector of names to the `covarianceModel`, 
where all elements in a block are specified as a comma separated list.

For example, given the desire for a block matrix between CL/V
and a diagnal for Ka, it would be coded as such:

```
modelSpec <- list(
  # Define the block diagonal elements
  covarianceModel = c("CL,V")
)
```

### `covariateModel`

The covariateModel element specifies the relationship between the parameter
and covariates of interest

```
modelSpec <- list(
  # ...
  covariateModel = list(
    CL = c("WT0","SEX"),
    ka = c("SEX"),
    Vc = c("WT0")
  )
)
```

### `covariateModelValues`

```
modelSpec <- list(
  # ...
  covariateModelValues = list(
    CL = c(WT0=0.75, SEX=0.3),
    ka = c(SEX=0.5),
    Vc = c(WT0=1)
  )
)
```

### `COVestimate`

Set covariate parameters to be estimated (1) or not (0)

```
modelSpec <- list(
  # ...
  COVestimate = list(
    CL = c(WT0=0, SEX=1),
    ka = c(SEX=1),
    Vc = c(WT0=0)
  )
)
```

### `COVcentering`

Define custom centering / reference values for covariates.
For continuous covariates, the default is to take the median 
value. (*note*: this behavior will be under discussion. This
was provided as-is for the initial implementation, however Devin
has significant issues with this approach. Namely, the transferability
of model interpretation is weak as the centering value will be different
across datasets. Furthermore, given a desire to selectively include/exclude
outlier subjects, the across the analysis a different centering value
may apply. Better to center on a reasonable reference value, for example, 
a 70kg, 40yr old male). 

for categorical values the default is the lowest group number

```
modelSpec <- list(
  # ...
  COVcentering  = c(WT0=70, RACE=2)
)
```

## algorithm spec
For the  arguments, an example may look like

```
algorithmSpec <- list(
  METHOD              = "FOCEI",
  MAXEVAL             = 9999,
  SIGDIGITS           = 3,
  PRINT               = 1,
  M4                  = FALSE,
  SEED                = 123456,
```

or for Iterative 2 stage into SAEM
```
algorithmSpec <- list(
  METHOD              = "SAEM",
  MAXEVAL             = 9999,
  SIGDIGITS           = 3,
  PRINT               = 1,
  M4                  = FALSE,
  SEED                = 123456,
  K1                  = 500,
  K2                  = 200,
  NRCHAINS            = 1,
  ITS                 = TRUE,
  ITS_ITERATIONS      = 50,
  IMPORTANCESAMPLING  = TRUE,
  IMP_ITERATIONS      = 50
)
```


## Additional Assumptions

- MU Referencing always used.
- Continuous covariates are always log-transformed, independent
  of the distribution of the parameter on which they are added. Centering
  by the median (or user defined values) for the covariates.
- Always untransformed categorical covariates. Several categories per
  covariate possible but all need to be numeric and integers.
- IIV correlation parameters always 0.1 at the initial guess.
- Selection of PRED,RES,WRES outputs dependent on the method that is used
  for estimation (in the tables renamed to: XPRED, XRES, XWRES):
     - PREDI RESI WRESI if FO
     - CPREDI, CRESI, CWRESI if FOCE
     - EPRED, ERES, EWRES if SAEM
- Default values for add and prop errors: 1 and 0.3
- Dataset can contain CMT or (ADM+YTYPE) columns. The output on the
  screen of this function will guide the user as to the needed values in
  these columns.
     - ADM+YTYPE but not CMT
         - YTYPE defines number of output
         - ADM used as CMT column (states reordered to accomodate)
     - ADM+CMT but not YTYPE
         - YTYPE is inferred based on CMT for observation records (in $ERROR)
           But this means that CMT needs to follow the OUTPUTn numbering!
         - CMT will be used as defined for selecting the dosing compartments
         - ADM is used to inform potential switchings for NONMEM parameters in the PK section


