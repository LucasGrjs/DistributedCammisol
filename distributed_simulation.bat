@echo off
TITLE CAMMISOL - DISTRIBUE (MPI)
SETLOCAL Enabledelayedexpansion

echo ==============================================================
echo    CONFIGURATION DE LA SIMULATION CAMMISOL (distribuee / MPI)
echo ==============================================================
echo.

:: 1. Lister les fichiers XML disponibles
echo Fichiers detectes dans Cammisol\models\cammisol\ :
echo --------------------------------------------------
if exist "Cammisol\models\cammisol\" (
    dir /b "Cammisol\models\cammisol\*.xml"
) else (
    echo [ERREUR] Dossier Cammisol\models\cammisol\ introuvable.
    pause
    exit /b
)
echo --------------------------------------------------
echo.

:: 2. Demander le nom du fichier
set /p NOM_FICHIER="Entrez le nom du fichier XML (ex: Distribution_CAMMISOL.xml) : "

:: Verification de l'existence du fichier
if not exist "Cammisol\models\cammisol\%NOM_FICHIER%" (
    echo.
    echo [ERREUR] Le fichier "%NOM_FICHIER%" est introuvable.
    pause
    exit /b
)

:: 3. Choix du nombre de processus MPI
set /p CORES="Nombre de processus MPI (Par defaut: 4) : "
if "%CORES%"=="" set CORES=4

echo.
echo [INFO] Preparation du conteneur...
echo [INFO] Fichier    : cammisol/%NOM_FICHIER%
echo [INFO] Processus  : %CORES%
echo.

:: 4. Lancement Docker (demarre le conteneur s'il n'est pas deja la)
docker-compose up -d

:: 5. Execution de la simulation
:: FILE_PATH est le chemin relatif a Cammisol/models utilise par startDistributedCammisol
set "FILE_PATH=cammisol/%NOM_FICHIER%"
echo [INFO] Lancement de la simulation dans le conteneur...
docker exec -it distributedcammisol_container bash -c "cd models && ./startDistributedCammisol %FILE_PATH% %CORES%"

echo.
echo Simulation terminee. Resultats dans Cammisol\models\output.log\
pause

docker-compose stop
