import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.io as pio
import base64
from pathlib import Path
from PIL import Image

pio.templates.default = "plotly_white"

# TODO:
# add number of detected DME´s a
# add Download CSV for detected DME´s
# add option to change start time for first BP


@st.cache
def load_data(path: str = "src/DB_MTX_USA.xlsx"):
    return pd.read_excel(path, parse_dates=["Sample_time"]).sort_values(
        ["Patient_id", "Sample_time"]
    )


# def load_data():
#     excel_file = st.sidebar.file_uploader("Choose your samples file", type=["xlsx"])

#     if excel_file is None:
#         st.info("Please specify samples and infusion time files in the sidebar")
#         return


def _to_number(x):
    """Remove the '<0.0x' strings from result..."""
    try:
        return float(x)
    except:
        return 0.0


def generate_download(df: pd.DataFrame) -> str:
    csv_file = df.to_csv(index=False, sep=";")
    b64 = base64.b64encode(csv_file.encode()).decode()
    return f'<a href="data:file/csv;base64,{b64}">Download CSV File</a> (right-click and save as &lt;some_name&gt;.csv)'


@st.cache
def read_markdown_file(markdown_file):
    return Path(markdown_file).read_text()


def main():

    image = Image.open("src/image/header.png")
    st.image(image, use_column_width=True)

    def display_sidebar_settings():
        desc_check = st.sidebar.checkbox("ℹ️ About")
        desc_markdown = read_markdown_file("src/desc_markdown.md")
        st.sidebar.markdown("---")
        st.sidebar.subheader("Configuration")

        if desc_check:
            st.sidebar.markdown(desc_markdown, unsafe_allow_html=True)

    display_sidebar_settings()

    # Load and sort data
    df = load_data()
    dose_mtx_df = df.loc[  # create a dose_mtx df
        df["Sample_type"] == "Dose_MTX", ["Patient_id", "Sample_time", "Result"]
    ]
    dose_mtx_df["Result"] = dose_mtx_df["Result"].apply(
        _to_number
    )  # apply to_number to remove the '<0.0x
    level_mtx_df = df.loc[  # create a level_mtx df
        df["Sample_type"] == "Level_MTX", ["Patient_id", "Sample_time", "Result"]
    ]
    level_mtx_df["Result"] = level_mtx_df["Result"].apply(
        _to_number
    )  # apply to_number to remove the '<0.0x

    # Set assumption for start treatment
    HOUR_FIRST_SAMPLE_TREATMENT = st.sidebar.slider(
        f"Assumption: first blood sample is ~hour:",
        0,
        48,
        23,
        1,
    )

    # Compute new treatment IDs
    start_treatment_threshold = st.sidebar.slider(
        "Start treatment threshold: ", 10, 150, 20, 5
    )
    level_mtx_df["next_result"] = level_mtx_df.groupby("Patient_id")["Result"].shift(-1)
    level_mtx_df["start_streak"] = (
        level_mtx_df["Result"] > start_treatment_threshold
    ) & (level_mtx_df["Result"] >= level_mtx_df["next_result"])
    level_mtx_df["treatment_id"] = (
        level_mtx_df.groupby("Patient_id")["start_streak"].cumsum().astype(str)
    )

    # Build hour difference within treatment
    level_mtx_df["hour_difference"] = level_mtx_df.groupby(
        ["Patient_id", "treatment_id"]
    )["Sample_time"].apply(lambda x: x - x.min()) / np.timedelta64(1, "h")
    level_mtx_df["hour_difference"] = (
        level_mtx_df["hour_difference"] + HOUR_FIRST_SAMPLE_TREATMENT
    )

    # Compute DME diagnostic
    hour_confidence_bound = st.sidebar.slider(
        "Select 42 + x hour for DME analysis: ", 1, 48, 4, 1
    )
    dme_hour_42_threshold = st.sidebar.slider(
        f"Conc. µM for DME between hour 42 -> 42+{hour_confidence_bound}: ",
        0.0,
        2.0,
        0.2,
        0.1,
    )
    level_mtx_df["is_dme"] = (
        (level_mtx_df["hour_difference"] >= 42)
        & (level_mtx_df["hour_difference"] < 42 + hour_confidence_bound)
        & (level_mtx_df["Result"] >= dme_hour_42_threshold)
    )
    is_patient_streak_dme = (
        level_mtx_df.groupby(["Patient_id", "treatment_id"])["is_dme"]
        .max()
        .reset_index(name="dme_diagnostic")
    )
    level_mtx_df = level_mtx_df.merge(
        is_patient_streak_dme, on=["Patient_id", "treatment_id"], how="left"
    )

    with st.beta_expander("Data preview"):
        st.dataframe(level_mtx_df.head(50))

    with st.beta_expander("Study particular patient"):
        all_patients = (
            level_mtx_df["Patient_id"]
            .value_counts()
            .sort_values(ascending=False)
            .index.values
        )
        number_unique_patients = st.write(
            f"Number of patients in data: {len(all_patients)}"
        )
        selected_patient = st.selectbox("Choose patient: ", all_patients)
        patient_level_df = level_mtx_df[level_mtx_df["Patient_id"] == selected_patient]
        patient_dose_df = dose_mtx_df[dose_mtx_df["Patient_id"] == selected_patient]
        patient_level_df["is_dose"] = False
        patient_dose_df["is_dose"] = True
        patient_df = patient_level_df.copy().append(patient_dose_df)

        use_log_x_scale = st.checkbox("Use logarithmic scale for y", True)
        fig_timeline = fig_per_sequence = px.scatter(
            patient_df,
            x="Sample_time",
            y="Result",
            color="is_dose",
            log_y=use_log_x_scale,
            title="Timeline",
        )
        st.plotly_chart(fig_timeline, use_container_width=True)
        use_hour_difference = st.checkbox(
            "Superimpose treatments per hour difference", False
        )
        fig_per_sequence = px.scatter(
            patient_level_df,
            x="hour_difference" if use_hour_difference else "Sample_time",
            y="Result",
            color="treatment_id",
            log_y=use_log_x_scale,
            title="Level per sequence",
        )
        st.plotly_chart(fig_per_sequence, use_container_width=True)

    with st.beta_expander("Study DME patients"):
        dme_patients = (
            level_mtx_df.loc[level_mtx_df["dme_diagnostic"] == True, "Patient_id"]
            .value_counts()
            .sort_values(ascending=False)
            .index.values
        )
        if len(dme_patients) == 0:
            st.info("No DME detected")
            return
        st.info(  # FIXME add number of DME detected
            f"Number of DME patients detected: {len(dme_patients)}"
        )
        selected_patient = st.selectbox("Choose patient: ", dme_patients)
        patient_level_df = level_mtx_df[level_mtx_df["Patient_id"] == selected_patient]

        if patient_level_df["treatment_id"].nunique() == 1:
            st.warning(
                "Beware, this patient only has one long treatment, maybe it's actually a mix of treatments"
            )
        st.plotly_chart(
            px.scatter(
                patient_level_df,
                x="hour_difference",
                y="Result",
                color="treatment_id",
                symbol="dme_diagnostic",
                title="Has patient DME?",
            ),
            use_container_width=True,
        )
        # FIXME: add why DME detection was triggered
        dme_patients_series = pd.Series(dme_patients)
        st.write(dme_patients_series)

    with st.beta_expander("Generate download link"):
        st.markdown(generate_download(dme_patients_series), unsafe_allow_html=True)


if __name__ == "__main__":
    st.set_page_config(
        page_title="DME Exploration",
        page_icon="chart",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    main()
