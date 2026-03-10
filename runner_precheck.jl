using YAML
using PackageScanner

# Dropbox Files On-Demand helper - force download entire directory
function force_download_directory(dirpath)
    @info "Forcing download of all files in directory (this triggers Dropbox sync)..." dirpath
    
    file_count = 0
    download_count = 0
    
    for (root, dirs, files) in walkdir(dirpath)
        for file in files
            filepath = joinpath(root, file)
            file_count += 1
            
            # Check if file is a stub (appears as 0 bytes)
            initial_size = filesize(filepath)
            
            # Force download by reading the entire file
            try
                # Reading the file will trigger Dropbox to download it
                open(filepath, "r") do io
                    # Read in chunks to avoid memory issues with large files
                    while !eof(io)
                        read(io, min(1024*1024, bytesavailable(io)))
                    end
                end
                
                # Check size after reading
                final_size = filesize(filepath)
                
                if initial_size == 0 && final_size > 0
                    download_count += 1
                    if download_count % 10 == 0
                        @info "Downloaded $download_count files so far..."
                    end
                end
            catch e
                @warn "Could not read file" filepath exception=e
            end
        end
    end
    
    @info "Processed $file_count files ($download_count were downloaded from Dropbox)"
    
    # Verify directory now has content
    size_output = chomp(read(`du -sh $dirpath`, String))
    @info "Final directory size: $size_output"
end

# Read configuration
vars = YAML.load_file(joinpath(ENV["GITHUB_WORKSPACE"], "_variables.yml"))

@info "Configuration loaded" vars

# Construct paths

source_path = joinpath(ENV["JPE_DBOX_APPS"], vars["dropbox_rel_path"], "replication-package")
dest_path = joinpath(ENV["GITHUB_WORKSPACE"], "replication-package")

@info "Paths configured" GITHUB_WORKSPACE=ENV["GITHUB_WORKSPACE"] source_path dest_path

# Check if source exists
@info "Checking source path exists..."
if !isdir(source_path)
    error("✗ Package not found at $source_path")
end
@info "✓ Source path exists"

# Check initial size
initial_size = chomp(read(`du -sh $source_path`, String))
@info "Initial package size (may be 0B if files are placeholders): $initial_size"

# Force download all files from Dropbox
@info "Ensuring all files are downloaded from Dropbox (this may take several minutes)..."
try
    force_download_directory(source_path)
catch e
    @error "Failed to download files" exception=e
    rethrow(e)
end

# Copy package
@info "Copying package to workspace..."
start_time = time()

try
    # Remove destination if it exists
    if isdir(dest_path)
        rm(dest_path; recursive=true, force=true)
    end
    @warn "using PackageScanner.mycp as workaround here"
    PackageScanner.mycp(source_path, dest_path; recursive = true, force=true)
    
    elapsed = time() - start_time
    @info "✓ Package copied successfully in $(round(elapsed, digits=2)) seconds"
catch e
    @error "Failed to copy package" exception=e
    rethrow(e)
end

# Unzip files depending on size
pkg_size = vars["package_size_gb"]
max_pkg_size = vars["package_max_pkg_size_gb"]
max_file_size = vars["package_max_file_size_gb"]
if pkg_size > max_pkg_size
    @info "package is larger than $(max_pkg_size) GB. Go into partial extraction mode"
    pkg_dir, manifest = PackageScanner.prepare_package_for_precheck(dest_path, size_threshold_gb = max_file_size, interactive = false)
    @info "Running precheck on $pkg_dir"

    PackageScanner.precheck_package(pkg_dir, pre_manifest=manifest)
else
    @info "Unzipping files in $dest_path"
    try
        zips = PackageScanner.read_and_unzip_directory(dest_path)
        @info "Unzipped $(length(zips)) file(s)"
    catch e
        @warn "Unzip had issues (may be okay)" exception=e
    end

    # Run precheck
    @info "Running precheck on $dest_path"
    try
        PackageScanner.precheck_package(dest_path)
        @info "✓ Precheck complete"
    catch e
        @error "Precheck failed" exception=e
        rethrow(e)
    end
end

