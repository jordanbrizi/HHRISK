::[Bat To Exe Converter]
::
::YAwzoRdxOk+EWAjk
::fBw5plQjdCyDJGyX8VAjFDFbXgCDMleeA6YX/Ofr08KvkGQuf6IpfZ3I07uAYLJdunnHWrII6lNXmcwJHydefwSibQA652dBuQQ=
::YAwzuBVtJxjWCl3EqQJgSA==
::ZR4luwNxJguZRRnk
::Yhs/ulQjdF+5
::cxAkpRVqdFKZSzk=
::cBs/ulQjdF+5
::ZR41oxFsdFKZSDk=
::eBoioBt6dFKZSDk=
::cRo6pxp7LAbNWATEpCI=
::egkzugNsPRvcWATEpCI=
::dAsiuh18IRvcCxnZtBJQ
::cRYluBh/LU+EWAnk
::YxY4rhs+aU+JeA==
::cxY6rQJ7JhzQF1fEqQJQ
::ZQ05rAF9IBncCkqN+0xwdVs0
::ZQ05rAF9IAHYFVzEqQJQ
::eg0/rx1wNQPfEVWB+kM9LVsJDGQ=
::fBEirQZwNQPfEVWB+kM9LVsJDGQ=
::cRolqwZ3JBvQF1fEqQJQ
::dhA7uBVwLU+EWDk=
::YQ03rBFzNR3SWATElA==
::dhAmsQZ3MwfNWATElA==
::ZQ0/vhVqMQ3MEVWAtB9wSA==
::Zg8zqx1/OA3MEVWAtB9wSA==
::dhA7pRFwIByZRRnk
::Zh4grVQjdCyDJGyX8VAjFDFbXgCDMleeA6YX/Ofr08KvkGQuf6IpfZ3I07uAYLJdunnHWrII6lJyquI/J1VdZhfL
::YB416Ek+ZG8=
::
::
::978f952a14a936cc963da21a135fa983
@echo off

color 0b

title HHRISK version 2

:main
cls
bg cursor 0
mode 110,30  

echo:
echo:
echo:
echo ************************************************************************************************************ 
echo ************************************************************************************************************
echo ****                                                                                                    ****
echo ****                                                                                                    ****
echo ****                 888    888  888    888  8888888b.   8888888  .d8888b.   888    d8P                 ****
echo ****                 888    888  888    888  888   Y88b    888   d88P  Y88b  888   d8P                  ****
echo ****                 888    888  888    888  888    888    888   Y88b.       888  d8P                   ****
echo ****                 8888888888  8888888888  888   d88P    888    "Y888b.    888d88K                    ****
echo ****                 888    888  888    888  8888888P"     888       "Y88b.  8888888b                   ****
echo ****                 888    888  888    888  888 T88b      888         "888  888  Y88b                  ****
echo ****                 888    888  888    888  888  T88b     888   Y88b  d88P  888   Y88b                 ****
echo ****                 888    888  888    888  888   T88b  8888888  "Y8888P"   888    Y88b                ****
echo ****                                                                                                    **** 
echo ****                                                                                                    ****
echo ************************************************************************************************************
echo ************************************************************************************************************  
echo:
echo:
echo:
echo:
echo                                                      MENU                                                                                                                                                           

insertbmp /p:"imagem.bmp" /x:0 /y:0 /z:30

cmdMenuSel f874 "How use HHRISK" "HHRISK Guide" "Input Data" "Execute HHRISK" "Quit"

if %ERRORLEVEL% == 1 goto starlink
if %ERRORLEVEL% == 2 goto starguide
if %ERRORLEVEL% == 3 goto starexcel
if %ERRORLEVEL% == 4 goto starcode
if %ERRORLEVEL% == 5 EXIT

:starexcel
START input.xlsm
goto main

:starcode
START RUNCODE
goto main

:starlink
START https://youtu.be/s2EJ9eF3Sj0
goto main

:starguide
START HHRISK_9_Guide.pdf
goto main
