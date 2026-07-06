# DistributedCammisol

# Distribution Model

This document outlines the **DistributionModel** project, which contains proof-of-concept implementations and experiments for distribution models within the **GAMA Platform**.

## Requirements

To set up and run this project, you'll only need:

* **Docker** and **Docker Compose**

Everything else (Java 17, Maven 3.8.6, Open MPI 4.1.4 with Java bindings, and the compiled GAMA/Cammisol build) is installed and built inside the Docker image, see [Dockerfile](Dockerfile).

## Building and Starting the Container

Follow these steps to build the image and start the container:

1.  Build the image (this compiles GAMA and Cammisol, so it might take a while):
    ```bash
    docker compose build
    ```
2.  Start the container in the background:
    ```bash
    docker compose up -d
    ```
3.  Open a shell inside the running container:
    ```bash
    docker compose exec app bash
    ```

All commands below are run from inside that shell, in the `/app/Cammisol` working directory.

---

## Running the CAMMISOL Model (for testing)

To start the CAMMISOL Model for testing purposes , execute the following command:



```bash
./startHeadless cammisol/cammisol.xml
 ```

## How to Start the Distribution Models

```bash
./startDistributedCammisol cammisol/Distribution_CAMMISOL.xml N # execute CAMMISOL on N process
```

---
### Results

All results from these simulations will be located in the `Cammisol/models/output.log/` directory after the model execution, it will contain the logs of the simulation from each **Processor**.

Since `./Cammisol` is bind-mounted into the container, the `models/` directory (including `.gaml` files and `output.log/`) can be edited/inspected directly from the host without rebuilding the image.

### Stopping the Container

```bash
docker compose down
```

