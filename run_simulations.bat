@echo off
setlocal enabledelayedexpansion

rem Runs a centralized CAMMISOL simulation, then a distributed one on 6 cores,
rem inside the docker compose "app" container, and finally compares their CO2
rem results using Cammisol/models/scripts/compare_co2.py on the host.

set CORES=6

echo === Making sure the container is up ===
docker compose up -d
if errorlevel 1 goto :error

echo.
echo === Running centralized simulation (cammisol.xml) ===
docker compose exec -T app bash -c "cd /app/Cammisol/models && ./startHeadless cammisol/cammisol.xml"
if errorlevel 1 goto :error

echo.
echo === Running distributed simulation (Distribution_CAMMISOL.xml, %CORES% cores) ===
docker compose exec -T app bash -c "cd /app/Cammisol/models && ./startDistributedCammisol cammisol/Distribution_CAMMISOL.xml %CORES%"
if errorlevel 1 goto :error

echo.
echo === Comparing CO2 results ===
python "%~dp0Cammisol\models\scripts\compare_co2.py"

endlocal
exit /b 0

:error
echo.
echo A step failed, aborting.
endlocal
exit /b 1
