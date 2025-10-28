@echo off
REM cLoupeComplement launcher for Windows
REM Requirements: Docker Desktop must be installed and running

set IMAGE=cloupecomplement

echo.
echo ========================================
echo   cLoupeComplement - Windows Launcher
echo ========================================
echo.

REM Check if Docker is running
docker info >nul 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Docker is not running!
    echo Please start Docker Desktop and try again.
    echo.
    pause
    exit /b 1
)

REM Check if image exists
docker image inspect %IMAGE% >nul 2>&1
if %errorlevel% neq 0 (
    echo Docker image '%IMAGE%' not found.
    echo.

    REM Check if pre-built image tar.gz exists
    if exist "docker\cloupecomplement-docker-image.tar.gz" (
        echo Found pre-built image file. Loading...
        echo This may take a few minutes...
        echo.
        gunzip -c docker\cloupecomplement-docker-image.tar.gz | docker load
        if %errorlevel% equ 0 (
            echo Image loaded successfully!
            echo.
        ) else (
            echo ERROR: Failed to load image from tar.gz
            echo The image file may be corrupted.
            echo.
            pause
            exit /b 1
        )
    ) else (
        echo ERROR: Docker image file not found!
        echo Please ensure 'docker\cloupecomplement-docker-image.tar.gz' exists.
        echo.
        pause
        exit /b 1
    )
)

echo Stopping any existing containers...
for /f "tokens=*" %%i in ('docker ps -q --filter ancestor=%IMAGE%') do docker stop %%i >nul 2>&1

echo Starting cLoupeComplement...
echo.

REM Start container in background
start /B docker run --rm ^
  -p 3838:3838 ^
  %IMAGE%

REM Wait for app to start
timeout /t 5 /nobreak >nul

REM Open browser
start http://localhost:3838

echo.
echo ========================================
echo   App Status
echo ========================================
echo.
echo ^>^> cLoupeComplement is running at: http://localhost:3838
echo ^>^> Upload your H5 and CSV files through the web interface
echo ^>^> Downloads will save to your browser's default download folder
echo.
echo IMPORTANT: Keep this window open!
echo Close this window to stop the app.
echo.
echo ========================================
echo.

REM Keep window open
pause
