create table if not exists sources (
    id text primary key,
    title text not null,
    revision text not null,
    url text not null
);

create table if not exists elements (
    atomic_number integer primary key,
    symbol text not null unique,
    name text not null,
    group_number integer,
    period_number integer not null,
    block text not null,
    category text not null,
    standard_state text,
    year_discovered text,
    atomic_weight text,
    atomic_weight_source text references sources(id)
);

create table if not exists compound_structures (
    id text primary key,
    canonical_formula text not null,
    charge integer not null,
    preferred_name text not null,
    canonical_smiles text not null,
    inchi text not null,
    inchi_key text not null,
    structure_kind text not null,
    source_id text not null references sources(id)
);
