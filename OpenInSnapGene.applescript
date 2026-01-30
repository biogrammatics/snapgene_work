-- OpenInSnapGene.applescript
-- Wrapper to open .dna files in SnapGene with a consistent window size
-- Save as Application in Script Editor, then use as default app for .dna files
-- or drag files/folders onto it

-- Configure your preferred window size here (left, top, right, bottom)
property windowBounds : {100, 100, 1400, 900}

on open theItems
    -- Collect all .dna files (expanding folders)
    set dnaFiles to {}
    repeat with anItem in theItems
        set itemInfo to info for anItem
        if folder of itemInfo then
            -- It's a folder - find all .dna files inside
            set folderPath to POSIX path of anItem
            set dnaFiles to dnaFiles & my findDNAFiles(folderPath)
        else
            -- It's a file - check if it's a .dna file
            if name of itemInfo ends with ".dna" then
                set end of dnaFiles to anItem
            end if
        end if
    end repeat

    if (count of dnaFiles) = 0 then
        display dialog "No .dna files found." buttons {"OK"} default button "OK"
        return
    end if

    tell application "SnapGene"
        activate
    end tell

    repeat with aFile in dnaFiles
        tell application "SnapGene"
            open aFile
        end tell

        -- Give SnapGene time to open the window
        delay 0.5

        -- Use System Events to resize the window
        tell application "System Events"
            tell process "SnapGene"
                if (count of windows) > 0 then
                    set position of front window to {item 1 of windowBounds, item 2 of windowBounds}
                    set size of front window to {(item 3 of windowBounds) - (item 1 of windowBounds), (item 4 of windowBounds) - (item 2 of windowBounds)}
                end if
            end tell
        end tell
    end repeat
end open

-- Find all .dna files in a folder (non-recursive for safety)
on findDNAFiles(folderPath)
    set dnaFiles to {}
    try
        -- Use find with pipe delimiter (do shell script strips newlines)
        set shellCmd to "find " & quoted form of folderPath & " -maxdepth 1 -name '*.dna' -type f 2>/dev/null | tr '\\n' '|'"
        set fileList to do shell script shellCmd
        if fileList is not "" then
            set oldDelims to AppleScript's text item delimiters
            set AppleScript's text item delimiters to "|"
            set fileLines to text items of fileList
            set AppleScript's text item delimiters to oldDelims
            repeat with filePath in fileLines
                if filePath is not "" then
                    set end of dnaFiles to (POSIX file filePath) as alias
                end if
            end repeat
        end if
    end try
    return dnaFiles
end findDNAFiles

-- Handle double-click on the app itself (no files)
on run
    display dialog "Drag .dna files onto this application to open them in SnapGene with a consistent window size." & return & return & "Current window size: " & (item 3 of windowBounds) - (item 1 of windowBounds) & " x " & ((item 4 of windowBounds) - (item 2 of windowBounds)) & " pixels" buttons {"OK"} default button "OK"
end run
