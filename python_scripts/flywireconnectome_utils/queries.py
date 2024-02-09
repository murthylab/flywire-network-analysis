from urllib.parse import urlparse
from sqlalchemy import MetaData, create_engine
from sqlalchemy.engine.url import make_url
from sqlalchemy.orm import scoped_session, sessionmaker
import os
import caveclient
import pandas as pd

HOME = os.path.expanduser("~")

datastack_name = "flywire_fafb_production"
cache_dir = f"{HOME}/FlyWireConnectome/query_cache"

if not os.path.exists(cache_dir):
    os.makedirs(cache_dir)


def set_env():
    os.environ[
        "SQLALCHEMY_DATABASE_URI"
    ] = "postgresql://postgresql:welcometothematrix@127.0.0.1:3306/annotation"
    os.environ["POSTGRES_USER"] = "postgres"
    os.environ["POSTGRES_PASSWORD"] = "welcometothematrix"
    os.environ["POSTGRES_DB"] = "annotation"
    os.environ["AUTH_URI"] = "global.daf-apis.com/auth"
    os.environ["AUTH_URL"] = "global.daf-apis.com/auth"
    os.environ["STICKY_AUTH_URL"] = "global.daf-apis.com/sticky_auth"
    os.environ["INFO_URL"] = "global.daf-apis.com/info"
    os.environ["DAF_CREDENTIALS"] = "/home/svenmd/.cloudvolume/secrets/cave-secret.json"


def create_session(sql_uri: str = None):
    engine = create_engine(
        sql_uri, pool_recycle=3600, pool_size=20, max_overflow=50, pool_pre_ping=True
    )
    Session = scoped_session(
        sessionmaker(bind=engine, autocommit=False, autoflush=False)
    )
    session = Session()
    return session, engine


def get_sql_url_params(sql_url):
    if not isinstance(sql_url, str):
        sql_url = str(sql_url)
    result = urlparse(sql_url)
    url_mapping = {
        "user": result.username,
        "password": result.password,
        "dbname": result.path[1:],
        "host": result.hostname,
        "port": result.port,
    }
    return url_mapping


def reflect_tables(sql_base, database_name):
    sql_uri = f"{sql_base}/{database_name}"
    engine = create_engine(sql_uri)
    meta = MetaData(engine)
    meta.reflect(views=True)
    tables = [table for table in meta.tables]
    engine.dispose()
    return tables


def ping_connection(session):
    is_database_working = True
    try:
        # to check database we will execute raw query
        session.execute("SELECT 1")
    except Exception as e:
        is_database_working = False
    return is_database_working


class SqlAlchemyCache:
    def __init__(self):
        self._engines = {}
        self._sessions = {}

    def get_engine(self, aligned_volume: str):
        if aligned_volume not in self._engines:
            SQL_URI_CONFIG = os.environ["SQLALCHEMY_DATABASE_URI"]
            sql_base_uri = SQL_URI_CONFIG.rpartition("/")[0]
            sql_uri = make_url(f"{sql_base_uri}/{aligned_volume}")
            self._engines[aligned_volume] = create_engine(
                sql_uri,
                pool_recycle=3600,
                pool_size=20,
                max_overflow=50,
                pool_pre_ping=True,
            )
        return self._engines[aligned_volume]

    def get(self, aligned_volume: str):
        if aligned_volume not in self._sessions:
            session = self._create_session(aligned_volume)
        session = self._sessions[aligned_volume]
        connection_ok = ping_connection(session)

        if not connection_ok:
            return self._create_session(aligned_volume)
        return session

    def _create_session(self, aligned_volume: str):
        engine = self.get_engine(aligned_volume)
        Session = scoped_session(sessionmaker(bind=engine))
        self._sessions[aligned_volume] = Session
        return self._sessions[aligned_volume]

    def invalidate_cache(self):
        self._engines = {}
        self._sessions = {}


def execute_query(query_str, mat_version, **kwargs):
    client = caveclient.CAVEclient(datastack_name, auth_token="7d795366041ad0cdcb6b42e9175c81d6")
    aligned_volume_name = client.info.get_datastack_info(datastack_name)[
        "aligned_volume"
    ]["name"]
    session = sqlalchemy_cache.get(aligned_volume_name)

    db_name = "{}__mat{}".format(datastack_name, mat_version)
    Session = sqlalchemy_cache.get(db_name)
    engine = sqlalchemy_cache.get_engine(db_name)

    query_df = pd.read_sql_query(query_str, con=engine.raw_connection(), **kwargs)

    return query_df


def per_neuron_neuropil_count(mat_version, synapse_side="pre", load_nt=False):
    set_env()

    if mat_version < 622:
        neuropil_version = 1
    else:
        neuropil_version = 2

    if synapse_side == "pre":
        syn_col = "pre_pt_root_id"
    else:
        syn_col = "post_pt_root_id"

    if neuropil_version == 2:
        neuropil_name = "fly_synapses_neuropil_v2"
        cache_path_wont = f"{cache_dir}/per_neuron_neuropilv2_filtered_count_{synapse_side}_{mat_version}.feather"
        cache_path_wnt = f"{cache_dir}/per_neuron_neuropilv2_filtered_count_ntavg_{synapse_side}_{mat_version}.feather"
    elif neuropil_version == 1:
        neuropil_name = "fly_synapses_neuropil"
        cache_path_wont = f"{cache_dir}/per_neuron_neuropil_filtered_count_{synapse_side}_{mat_version}.feather"
        cache_path_wnt = f"{cache_dir}/per_neuron_neuropil_filtered_count_ntavg_{synapse_side}_{mat_version}.feather"
    else:
        raise Exception("Neuropil version unknown")

    if load_nt:
        cache_path = cache_path_wnt
    else:
        cache_path = cache_path_wont

    if os.path.exists(cache_path.replace("filtered_", "")):
        return pd.read_feather(cache_path.replace("filtered_", ""))

    nts = ["gaba", "ach", "glut", "oct", "ser", "da"]
    nt_str = ",".join(f"AVG({nt}) as {nt}_avg" for nt in nts)
    nt_cols = [f"{nt}_avg" for nt in nts]

    if not os.path.exists(cache_path):
        if load_nt:
            query_str = f"SELECT {syn_col}, neuropil, COUNT(*), {nt_str} FROM synapses_nt_v1 INNER JOIN valid_synapses_nt_v2 ON synapses_nt_v1.id=valid_synapses_nt_v2.target_id LEFT JOIN fly_synapses_neuropil_v2 ON synapses_nt_v1.id=fly_synapses_neuropil_v2.target_id WHERE pre_pt_root_id != 0 AND post_pt_root_id != 0 AND pre_pt_root_id != post_pt_root_id AND cleft_score > 50 GROUP BY ({syn_col}, neuropil)"
            # query_str = f"SELECT {syn_col}, neuropil, COUNT(*), {nt_str} FROM synapses_nt_v1 JOIN fly_synapses_neuropil ON synapses_nt_v1.id=fly_synapses_neuropil.target_id WHERE pre_pt_root_id != 0 AND post_pt_root_id != 0 AND pre_pt_root_id != post_pt_root_id AND cleft_score > 50 GROUP BY ({syn_col}, neuropil)"
            query_df = execute_query(query_str, mat_version)

            query_df.drop(columns=nt_cols).to_feather(
                cache_path_wont, compression="zstd"
            )
            query_df.to_feather(cache_path_wnt, compression="zstd")
        else:
            query_str = f"SELECT {syn_col}, neuropil, COUNT(*) FROM synapses_nt_v1 INNER JOIN valid_synapses_nt_v2 ON synapses_nt_v1.id=valid_synapses_nt_v2.target_id LEFT JOIN fly_synapses_neuropil_v2 ON synapses_nt_v1.id=fly_synapses_neuropil_v2.target_id WHERE pre_pt_root_id != 0 AND post_pt_root_id != 0 AND pre_pt_root_id != post_pt_root_id AND cleft_score > 50 GROUP BY ({syn_col}, neuropil)"
            # query_str = f"SELECT {syn_col}, neuropil, COUNT(*) FROM synapses_nt_v1 JOIN fly_synapses_neuropil ON synapses_nt_v1.id=fly_synapses_neuropil.target_id WHERE pre_pt_root_id != 0 AND post_pt_root_id != 0 AND pre_pt_root_id != post_pt_root_id AND cleft_score > 50 GROUP BY ({syn_col}, neuropil)"
            query_df = execute_query(query_str, mat_version)

            query_df.to_feather(cache_path_wont, compression="zstd")

    query_df = pd.read_feather(cache_path)

    return query_df


sqlalchemy_cache = SqlAlchemyCache()
