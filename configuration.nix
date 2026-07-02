{ config, pkgs, lib, ... }:

let
  # Runtime libraries for unpatched/generic Linux binaries (binary Python
  # wheels, downloaded toolchains, vendor blobs) through nix-ld.
  # Keep this list reasonably broad: if Claude Code installs something that
  # needs a library not listed here, the failure mode is a cryptic
  # `cannot open shared object file`, and the fix requires editing this file
  # *and* a full rebuild — annoying mid-session.
  runtimeLibs = with pkgs; [
    stdenv.cc.cc.lib  # libstdc++.so.6, libgcc_s.so.1
    gfortran.cc.lib   # libgfortran.so.5
    glibc

    blas
    lapack

    zlib
    zstd
    bzip2
    xz
    openssl
    libffi
    sqlite
    expat

    curl
    ncurses
    util-linux
  ];
in
{
  # Allow Claude Code (and any unfree packages it might want to add).
  nixpkgs.config.allowUnfree = true;

  # For easier rebuilding
  nix.channel.enable = true;
  nix.settings.nix-path = [ "nixpkgs=channel:nixos-unstable" ];

  # Point /etc/nixos/configuration.nix at the copy in the shared repo so
  # `sudo nixos-rebuild` picks up edits made on the host without any manual
  # symlinking. The target is a literal path string, so the file is read
  # live from the shared mount (not copied into the Nix store).
  environment.etc."nixos/configuration.nix".source =
    "/mnt/host/two-heavy-atoms-electronegativity/configuration.nix";

  # On a freshly built VM the root channel is registered but the tarball
  # hasn't been fetched, so `nixos-rebuild` fails with "file 'nixpkgs/nixos'
  # was not found". Fetch it on first activation if missing.
  system.activationScripts.bootstrapNixChannel = ''
    if [ ! -e /nix/var/nix/profiles/per-user/root/channels/nixos ]; then
      ${pkgs.nix}/bin/nix-channel --update || true
    fi
  '';

  boot.loader.grub.enable = false;

  virtualisation.graphics = false;

  # https://github.com/NixOS/nixpkgs/issues/499166
  documentation.doc.enable = false;

  # Vendor binaries with hardcoded shebangs (#!/bin/bash, etc.) resolve via
  # envfs. Without this, lots of `npm install`'d CLIs, vendor SDKs, and
  # ad-hoc scripts will fail on NixOS.
  services.envfs.enable = true;

  virtualisation.diskSize = 64 * 1024; # MiB — generous so `nix-store` growth
                                       # mid-session doesn't wedge a rebuild.
  # Default puts the writable store overlay on tmpfs, which steals from the
  # VM's RAM. With only 8 GiB, large rebuilds (quarto + R + texlive) get
  # OOM-killed. Put the overlay on the VM disk instead.
  virtualisation.writableStoreUseTmpfs = false;
  virtualisation.memorySize = 8 * 1024; # MiB
  virtualisation.cores = 4;

  virtualisation.sharedDirectories.vmShare = {
    source = "/home/alexlm/vm-shared";
    target = "/mnt/host";
    securityModel = "none";
  };

  imports = [
    <nixpkgs/nixos/modules/virtualisation/qemu-vm.nix>
  ];

  # SSH server, key-only.
  services.openssh = {
    enable = true;
    settings = {
      PasswordAuthentication = false;
      KbdInteractiveAuthentication = false;
      PermitRootLogin = "no";
    };
  };

  # LAN-accessible:
  #   ssh -i ~/.ssh/id_alexlm -p 2222 alexlm@<host-LAN-IP>
  virtualisation.forwardPorts = [
    {
      from = "host";
      host.address = "0.0.0.0";
      host.port = 2222;
      guest.port = 22;
    }
  ];

  users.users.alexlm = {
    isNormalUser = true;
    description = "Alex";
    extraGroups = [ "wheel" ];
    openssh.authorizedKeys.keys = [
      "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIN/c+S1Osh+QIBAQBy3iWB3n134WrFcZhdzOROJuf29k alexlovesmolecules@gmail.com"
      "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIP7faEDi8p7cjUMNNfx3Ir/B6FDAPgXzaLod35j/LI3/ alexlm@Alexanders-MacBook-Air.local"
    ];
  };

  # Claude Code needs to run `sudo nixos-rebuild switch` without an
  # interactive password prompt — agents can't type passwords.
  security.sudo.wheelNeedsPassword = false;

  # nix-ld for binary wheels and downloaded toolchains.
  programs.nix-ld = {
    enable = true;
    libraries = runtimeLibs;
  };

  # Make flakes + the new CLI available. Claude Code edits & rebuilds work
  # without flakes too, but they're the norm in modern docs/examples, so
  # enabling them avoids steering Claude into the legacy channel path
  # mid-session.
  nix.settings.experimental-features = [ "nix-command" "flakes" ];

  # Don't let the store fill up the disk after a few rebuilds.
  nix.gc = {
    automatic = true;
    dates = "weekly";
    options = "--delete-older-than 14d";
  };
  nix.settings.auto-optimise-store = true;

  environment.sessionVariables = {
    LD_LIBRARY_PATH = lib.makeLibraryPath runtimeLibs;
    LIBRARY_PATH = lib.makeLibraryPath (with pkgs; [ blas lapack ]);

    # If Claude installs Python tooling via uv, this keeps it on the nixpkgs
    # interpreter rather than downloading generic builds.
    UV_PYTHON_DOWNLOADS = "never";
    UV_PYTHON_PREFERENCE = "only-system";
  };

  environment.systemPackages = with pkgs; [
    bashInteractive  # envfs needs bash on PATH

    # The actual point of this config.
    claude-code

    # Editing + rebuilding the system. `nixos-rebuild` is in the system
    # closure automatically, so it's not listed here.
    vim
    git

    # General-purpose tools Claude reaches for constantly. Putting them in
    # the base image avoids a rebuild-roundtrip the first time Claude wants
    # to grep a codebase.
    ripgrep
    fd
    jq
    curl
    wget
    unzip
    file
    tree
    which

    # Build toolchain — pre-staged so source builds (`uv sync`, `npm install`,
    # `cargo build`, etc.) work without a rebuild.
    pkg-config
    gcc
    gnumake
    gfortran

    # Language runtimes Claude is most likely to want. Each one swapped in
    # via `nixos-rebuild` would otherwise be a 30-60 second interruption.
    python3
    uv
    nodejs
    cargo
    rustc
    go
  ];

  system.stateVersion = "25.11";
}
