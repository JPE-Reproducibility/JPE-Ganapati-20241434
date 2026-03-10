close all

e_bar   = 1.3
alpha_e = 3.0
alpha_c = 0.4

plot_theta_e_vs_log_n(alpha_e, e_bar);
saveas(gcf,'../../Output/F1_plot_theta_e_vs_log_n','epsc')

plot_theta_c_vs_log_n(alpha_c);
saveas(gcf,'../../Output/F1_plot_theta_c_vs_log_n','epsc')

plot_theta_combined(alpha_c, alpha_e, e_bar);
saveas(gcf,'../../Output/F1_plot_theta_vs_log_n','epsc')

plot_theta_i_vs_log_n(alpha_c, alpha_e, e_bar);
saveas(gcf,'../../Output/F1_plot_theta_i_vs_log_n','epsc')


generate_pdf_entry(e_bar, alpha_e);
saveas(gcf,'../../Output/AF1_pdf_e','epsc')

plot_Er_vs_log_n(alpha_c);
saveas(gcf,'../../Output/AF1_fEr','epsc')
