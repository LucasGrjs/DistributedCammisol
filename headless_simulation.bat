@echo off
TITLE CAMMISOL - CENTRALISE
SETLOCAL Enabledelayedexpansion

echo ==================================================
echo    CONFIGURATION DE LA SIMULATION CAMMISOL (centralisee)
echo ==================================================
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
set /p NOM_FICHIER="Entrez le nom du fichier XML (ex: cammisol.xml) : "

:: Verification de l'existence du fichier
if not exist "Cammisol\models\cammisol\%NOM_FICHIER%" (
    echo.
    echo [ERREUR] Le fichier "%NOM_FICHIER%" est introuvable.
    pause
    exit /b
)

echo.
echo [INFO] Preparation du conteneur...
echo [INFO] Fichier : cammisol/%NOM_FICHIER%
echo.

:: 3. Lancement Docker (demarre le conteneur s'il n'est pas deja la)
docker-compose up -d

:: 4. Execution de la simulation
:: FILE_PATH est le chemin relatif a Cammisol/models utilise par startHeadless
set "FILE_PATH=cammisol/%NOM_FICHIER%"
docker exec -it distributedcammisol_container bash -c "cd models && ./startHeadless %FILE_PATH%"

echo.
echo Simulation terminee. Resultats dans Cammisol\models\output.log\
pause

docker-compose stop
