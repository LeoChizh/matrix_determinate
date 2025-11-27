import os

print("Reading project files...")

file_count = 0
line_count = 0

with open("snapshot.txt", "w") as out:
    for root, dirs, files in os.walk("."):
        # Skip all unnecessary folders
        for skip in ['build', '.git', 'utils', 'fonts', 'textures', 'KHR', 'glad','.vscode']:
            if skip in dirs: 
                dirs.remove(skip)
        
        for file in files:
            # Skip specific files and binary files
            if file in ['glad.c', 'stb_image.h', 'snapshot.py', 'snapshot.txt', '.gitignore'] or file.endswith(('.exe', '.dll', '.obj', '.md')):
                continue
                
            path = os.path.join(root, file)
            file_count += 1
            
            try:
                if os.path.getsize(path) < 1000000:  # 1MB limit
                    with open(path, 'r') as f:
                        content = f.read()
                        file_lines = content.count('\n') + (1 if content else 0)
                        line_count += file_lines
                        print(f"Reading: {path} ({file_lines} lines)")  # Show file with line count
                        
                        out.write(f"\n=== {path} ===\n")
                        out.write(content)
                else:
                    print(f"Reading: {path} (skipped - file too large)")
            except: 
                print(f"Reading: {path} [ERROR reading file]")
                pass

print(f"\nDone! Processed {file_count} files with {line_count} total lines. Created snapshot.txt")