@echo off
setlocal enabledelayedexpansion

set "output=mdsp-ii.md"

:: 删除旧的合并文件，确保每次都是全新生成
if exist "%output%" del "%output%"

set "isFirst=1"
for /f "delims=" %%f in ('dir /b /on *.md') do (
    if /i not "%%f"=="%output%" if /i not "%%f"=="README.md" (
        if "!isFirst!"=="1" (
            type "%%f" >> "%output%"
            set "isFirst=0"
        ) else (
            :: 追加正文
            type "%%f" >> "%output%"
        )
    )
)

echo done %output%