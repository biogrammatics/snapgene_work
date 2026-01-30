-- OpenInSnapGene.applescript
-- Wrapper to open .dna files in SnapGene with a consistent window size
-- Save as Application in Script Editor, then use as default app for .dna files
-- or drag files onto it

-- Configure your preferred window size here (left, top, right, bottom)
property windowBounds : {100, 100, 1400, 900}

on open theFiles
    tell application "SnapGene"
        activate
    end tell

    repeat with aFile in theFiles
        set filePath to POSIX path of aFile

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

-- Handle double-click on the app itself (no files)
on run
    display dialog "Drag .dna files onto this application to open them in SnapGene with a consistent window size." & return & return & "Current window size: " & (item 3 of windowBounds) - (item 1 of windowBounds) & " x " & ((item 4 of windowBounds) - (item 2 of windowBounds)) & " pixels" buttons {"OK"} default button "OK"
end run
