
function download_and_check_md5sum(url::String, filename::String; md5_goal::String = "" )
    need_to_download = true
    if isfile(download_filename)
        if Sys.which("md5sum") == nothing
            println("# System doesn't have md5sum so forcing redownload of " * filename )
        else
            md5_here = ""
            try
                md5sum_output = read(`md5sum $filename`, String)
                md5_here = split(md5sum_output)[1]
                if length(md5_goal) == 0 || md5_here == md5_goal
                    println("# Found ", filename, " so not downloading.")
                    need_to_download = false
                else
                    println("# Found ", filename, " but md5sum doesn't match, so should redownload.")
                end
            catch
                println("# Warning couldn't compute md5sum to verify download of ", filename)
                need_to_download = false
            end
        end
    end
    if need_to_download
        println("# Downloading ", filename)
        download(url, filename)
        if Sys.which("md5sum") == nothing
            @warn "System doesn't have md5sum to check that download of " * filename * " was successful."
        else
            md5sum_output = read(`md5sum $filename`, String)
            md5_here = split(md5sum_output)[1]
            if md5_here != md5_goal
                println("# The md5sum didn't match expectations, make sure the download was successfuly completed.")
            end
        end
    end
    return
end
