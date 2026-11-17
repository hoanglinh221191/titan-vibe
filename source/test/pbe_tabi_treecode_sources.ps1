function Get-PbeTabiTreecodeSources {
    param([Parameter(Mandatory = $true)][string]$RepoRoot)

    @(
        (Join-Path $RepoRoot 'src\tabi\tabi_treecode.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_treecode_internal.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_treecode_geometry.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_treecode_series.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_treecode_policy.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_treecode_leaf.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_treecode_operator.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_fmm_m2l_prototype.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_fmm_m2l_fallback.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_fmm_m2l_cartesian.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_fmm_m2l_local.cpp')
        # tabi_near.cpp, which every caller compiles, calls curved_node_pair.
        (Join-Path $RepoRoot 'src\tabi\tabi_node_patch_curved.cpp')
    )
}

# What src\tabi\tabi_surface.cpp links against beyond the list above: the
# Calderon/SES surface builders it calls and their own dependencies.
function Get-PbeTabiSurfaceSources {
    param([Parameter(Mandatory = $true)][string]$RepoRoot)

    @(
        (Join-Path $RepoRoot 'src\tabi\tabi_calderon_aux_surface.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_calderon_trace_space.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_analytic_ses.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_ses_bilevel_grid.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_intersections.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_provenance.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_ses_patch_ray.cpp')
        (Join-Path $RepoRoot 'src\tabi\tabi_ses_quartic.cpp')
    )
}
