# Building and Pushing EvolvDX Exome Container Images

This guide explains how to build, tag, and push the custom EvolvDX exome container images to Amazon ECR for use with the pipeline on AWS Batch.

## Prerequisites

- Docker installed and running
- AWS CLI configured with appropriate credentials
- Access to the EvolvDX ECR registry (account ID: 062689759574, region: us-east-1)

## Included Container Images

The pipeline uses the following custom container images (stored in `env/exo/*/`):
- `exosomedx/methylkit`: R environment for methylKit QC reporting
- `exosomedx/basic_read_statistics`: Custom read statistics tools
- `exosomedx/methsnsv`: SNV calling from methylation data
- `exosomedx/custom_multiqc_picardhs`: MultiQC with Picard HS metrics support

## Build Workflow

### 1. Build the Docker image
> [!IMPORTANT]
> AWS Batch uses x86_64 (amd64) architecture. If you are building on a Mac with Apple Silicon (arm64), you must specify the `--platform linux/amd64` flag to ensure compatibility.

```bash
# Navigate to the directory containing the Dockerfile
cd env/exo/<package_name>/

# Build the image (replace <package_name> and <version> with appropriate values)
docker build --no-cache --platform linux/amd64 -t exosomedx/<package_name>:<version> .
```

### 2. Verify the image was built
```bash
docker image ls | grep exosomedx/<package_name>
```

### 3. Tag the image for ECR
```bash
# Get the IMAGE ID from the previous step, or use the tag
docker tag exosomedx/<package_name>:<version> 062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/<package_name>:<version>
```

### 4. Authenticate with ECR
```bash
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 062689759574.dkr.ecr.us-east-1.amazonaws.com
```

### 5. Create the ECR repository (if it doesn't exist)
```bash
aws ecr create-repository --repository-name exosomedx/<package_name> --region us-east-1
```

### 6. Push the image to ECR
```bash
docker push 062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/<package_name>:<version>
```

### 7. Update the pipeline config
After pushing the new image version, update the container path in the corresponding module file (`modules/local/exo/<package_name>.nf`) to reference the new image tag:
```nextflow
container "062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/<package_name>:<version>"
```

## Example: Building the methylKit Image

```bash
# Build
cd env/exo/methylkit/
docker build --no-cache --platform linux/amd64 -t exosomedx/methylkit:1.36.0 .

# Tag
docker tag exosomedx/methylkit:1.36.0 062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.36.0

# Push
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 062689759574.dkr.ecr.us-east-1.amazonaws.com
docker push 062689759574.dkr.ecr.us-east-1.amazonaws.com/exosomedx/methylkit:1.36.0

# Update module config
sed -i '' 's|exosomedx/methylkit:.*|exosomedx/methylkit:1.36.0|' modules/local/exo/methylkit.nf
```

## Versioning Convention

Use semantic versioning for container images:
- Major version (X.0.0): Breaking changes to the container environment
- Minor version (x.Y.0): New tools or features added, backward compatible
- Patch version (x.y.Z): Bug fixes, tool updates, backward compatible

## Troubleshooting

### "No space left on device" during build
Clean up unused Docker images and containers:
```bash
docker system prune -a
```

### "Denied: Your authorization token has expired" when pushing
Re-authenticate with ECR:
```bash
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 062689759574.dkr.ecr.us-east-1.amazonaws.com
```

### Image works locally but fails on AWS Batch
Ensure you built the image with `--platform linux/amd64` if you are on an ARM-based Mac.
