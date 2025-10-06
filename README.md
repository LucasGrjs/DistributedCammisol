# DistributedCammisol

# Distribution Model

This document outlines the **DistributionModel** project, which contains proof-of-concept implementations and experiments for distribution models within the **GAMA Platform**.

## Requirements

To set up and run this project, you'll need the following:

* **Java 17**
* **mpirun (Open MPI) 4.1.4**: Download from [https://www.open-mpi.org/software/ompi/v4.1/](https://www.open-mpi.org/software/ompi/v4.1/)
* **Apache Maven 3.8.6**
* **Java Binding for Open MPI**: Refer to [https://www.open-mpi.org/faq/?category=java](https://www.open-mpi.org/faq/?category=java) for details.

## Compiling the Project

Follow these steps to compile the project:

1.  Navigate to the `gama` directory:
    ```bash
    cd gama
    ```
2.  Run the build script. This process might take some time:
    ```bash
    ./travis/build.sh
    ```
3.  Change to the `Cammisol` directory:
    ```bash
    cd ../Cammisol/models
    ```
4.  Give execution permissions to the script:
    ```bash
    chmod +x startHeadless
    chmod +x startDistributedCammisol
    ```

---

## Running the CAMMISOL Model (for testing)

To start the CAMMISOL Model for testing purposes , execute the following command:



```bash
./startHeadless cammisol/cammisol.xml # execute the CAMMISOL Model on 1 process 
 ```
Warning : l.136 (do die;) of cammisol.gaml must be uncommented for centralized execution but commented for distributed execution

## How to Start the Distribution Models

```bash
./startDistributedCammisol cammisol/Distribution_CAMMISOL.xml.xml N # execute CAMMISOL on N process
```

---
### Results

All results from these simulations will be located in the `/output.log/` directory after the model execution, it will contain the logs of the simulation from each **Processor**.

