# Manual Push Instructions

Date: 2026-07-03

The public archive repository has been staged and committed locally.

Local repository:

```text
/home/lzb/glioma/Project_Management/Submission_Package/GitHub_Public_Archive_20260703
```

Remote repository:

```text
git@github.com:lzbkk/glioma-public-data-archive.git
https://github.com/lzbkk/glioma-public-data-archive
```

## Pre-Push Checks

Run:

```bash
cd /home/lzb/glioma/Project_Management/Submission_Package/GitHub_Public_Archive_20260703
git status
sha256sum -c CHECKSUMS.sha256
```

Expected:

```text
working tree clean
archive/Source_Data_Package_20260702.tar.gz: OK
archive/Journal_Supplementary_Tables_20260702.tar.gz: OK
archive/Supplementary_Tables_1_11_JournalFacing.xlsx: OK
```

## Push

If the GitHub repository is empty:

```bash
cd /home/lzb/glioma/Project_Management/Submission_Package/GitHub_Public_Archive_20260703
git push -u origin main
```

If GitHub was created with an initial README or other file and push is rejected,
do not force-push immediately. First inspect:

```bash
git fetch origin
git log --oneline --graph --decorate --all --max-count=20
```

Then decide whether to merge/rebase or recreate the remote repository as empty.

## Public Disclosure Boundary

Before pushing, confirm that public disclosure is intended. This archive excludes:

- raw GBM-Space H5AD/tarball files;
- raw GLASS Synapse files;
- raw third-party omics matrices;
- `.Renviron`, API keys, tokens, cookies and proxy settings;
- large local caches and private logs.

The archive includes derived source-data packages, supplementary workbooks,
analysis scripts and metadata/checksum files.

## After Push

After the GitHub push succeeds:

1. Create or connect the Zenodo record.
2. Confirm final licence choices.
3. Create a GitHub release if using the GitHub-Zenodo integration.
4. Record the Zenodo DOI and record URL.
5. Backfill the DOI into Data/Code Availability and project management files.
