# Results Page
elif page == "Results":
    st.title("Analysis Results")
    
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        
        # Display DataFrame
        st.dataframe(results_df)

        # Motif Occurrence Summary
        motif_counts = results_df["Motif"].value_counts().reset_index()
        motif_counts.columns = ["Motif", "Count"]
        st.subheader("Motif Occurrence Summary")
        st.dataframe(motif_counts)

        # Convert DataFrame to CSV
        csv = results_df.to_csv(index=False).encode('utf-8')

        # Download CSV button
        st.download_button(
            label="Download CSV",
            data=csv,
            file_name="motif_analysis_results.csv",
            mime="text/csv"
        )

    else:
        st.warning("No results available. Please upload or paste a sequence first.")
