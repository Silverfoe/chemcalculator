package com.chemcalc.ui;

import java.awt.BorderLayout;
import java.awt.Component;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.Objects;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.BiConsumer;
import javafx.application.Platform;
import javafx.embed.swing.SwingNode;
import javafx.scene.Node;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import me.friwi.jcefmaven.CefAppBuilder;
import org.cef.CefApp;
import org.cef.CefClient;
import org.cef.CefSettings;
import org.cef.browser.CefBrowser;
import org.cef.browser.CefFrame;
import org.cef.handler.CefDisplayHandlerAdapter;
import org.cef.handler.CefLoadHandlerAdapter;

final class ChromiumStructureViewer {
    private static final String BRIDGE_PREFIX = "__CHEMCALC__:";
    private static final Path INSTALL_DIR = Path.of(System.getProperty("user.home"), ".chemcalculator", "jcef-bundle");
    private static final Path ASSET_DIR = Path.of(System.getProperty("java.io.tmpdir"), "chemcalculator", "structure-viewer");

    private static CefApp cefApp;

    private final SwingNode swingNode = new SwingNode();
    private final JPanel hostPanel = new JPanel(new BorderLayout());
    private final BiConsumer<String, String> onAtomSelected;
    private final AtomicBoolean pageLoaded = new AtomicBoolean(false);

    private CefClient client;
    private CefBrowser browser;
    private String viewerUrl;
    private String pendingScript = StructureViewerHtmlRenderer.renderMessageScript(
            "No structure loaded",
            "Load a formula, SMILES, Molfile/SDF, or InChI to populate the viewer.");
    private String pendingStylePreset = "ball-stick";

    ChromiumStructureViewer(BiConsumer<String, String> onAtomSelected) {
        this.onAtomSelected = Objects.requireNonNull(onAtomSelected, "onAtomSelected");
        Platform.runLater(() -> swingNode.setContent(hostPanel));
        initialize();
    }

    Node node() {
        return swingNode;
    }

    void renderScript(String script, String stylePreset) {
        pendingScript = script;
        pendingStylePreset = stylePreset == null ? "ball-stick" : stylePreset;
        executePendingState();
    }

    void showFallbackMessage(String title, String detail) {
        pendingScript = StructureViewerHtmlRenderer.renderMessageScript(title, detail);
        executePendingState();
    }

    void setStylePreset(String stylePreset) {
        pendingStylePreset = stylePreset == null ? "ball-stick" : stylePreset;
        executeJavaScript("window.setStylePreset(" + quoteJs(pendingStylePreset) + ");");
    }

    void dispose() {
        SwingUtilities.invokeLater(() -> {
            if (browser != null) {
                browser.close(true);
            }
        });
    }

    private void initialize() {
        SwingUtilities.invokeLater(() -> {
            try {
                viewerUrl = prepareAssetDirectory().resolve("structure-viewer.html").toUri().toString();
                CefApp app = getOrCreateApp();
                client = app.createClient();
                client.addDisplayHandler(new CefDisplayHandlerAdapter() {
                    @Override
                    public boolean onConsoleMessage(
                            CefBrowser activeBrowser,
                            CefSettings.LogSeverity level,
                            String message,
                            String source,
                            int line) {
                        handleBridgeMessage(message);
                        return false;
                    }
                });
                client.addLoadHandler(new CefLoadHandlerAdapter() {
                    @Override
                    public void onLoadEnd(CefBrowser activeBrowser, CefFrame frame, int httpStatusCode) {
                        if (!frame.isMain()) {
                            return;
                        }
                        pageLoaded.set(true);
                        executePendingState();
                    }
                });

                browser = client.createBrowser(viewerUrl, false, false);
                Component component = browser.getUIComponent();
                hostPanel.removeAll();
                hostPanel.add(component, BorderLayout.CENTER);
                hostPanel.revalidate();
                hostPanel.repaint();
            } catch (Exception exception) {
                installFailurePane(exception);
            }
        });
    }

    private void handleBridgeMessage(String message) {
        if (message == null || !message.startsWith(BRIDGE_PREFIX)) {
            return;
        }
        String payload = message.substring(BRIDGE_PREFIX.length());
        if (payload.contains("\"type\":\"atom-selected\"")) {
            String atomIndex = extractJsonValue(payload, "atomIndex");
            String label = extractJsonValue(payload, "label");
            Platform.runLater(() -> onAtomSelected.accept(atomIndex, label));
        }
    }

    private void executePendingState() {
        executeJavaScript("window.setStylePreset(" + quoteJs(pendingStylePreset) + ");");
        executeJavaScript(pendingScript);
        executeJavaScript("window.setStylePreset(" + quoteJs(pendingStylePreset) + ");");
    }

    private void executeJavaScript(String script) {
        SwingUtilities.invokeLater(() -> {
            if (!pageLoaded.get() || browser == null) {
                return;
            }
            browser.executeJavaScript(script, viewerUrl, 0);
        });
    }

    private void installFailurePane(Exception exception) {
        hostPanel.removeAll();
        JLabel label = new JLabel(
                "<html><div style='padding:18px;font-family:Segoe UI;'>Chromium viewer failed to start.<br/>"
                        + escapeHtml(exception.getMessage())
                        + "</div></html>",
                SwingConstants.CENTER);
        hostPanel.add(label, BorderLayout.CENTER);
        hostPanel.revalidate();
        hostPanel.repaint();
    }

    private static synchronized CefApp getOrCreateApp() throws Exception {
        if (cefApp != null) {
            return cefApp;
        }
        Files.createDirectories(INSTALL_DIR);

        CefAppBuilder builder = new CefAppBuilder();
        builder.setInstallDir(INSTALL_DIR.toFile());
        builder.getCefSettings().windowless_rendering_enabled = false;
        builder.getCefSettings().background_color = builder.getCefSettings().new ColorType(255, 255, 255, 255);
        cefApp = builder.build();
        return cefApp;
    }

    private static Path prepareAssetDirectory() throws IOException {
        Files.createDirectories(ASSET_DIR.resolve("vendor"));
        copyResource("/ui/structure-viewer.html", ASSET_DIR.resolve("structure-viewer.html"));
        copyResource("/ui/structure-viewer.css", ASSET_DIR.resolve("structure-viewer.css"));
        copyResource("/ui/structure-viewer.js", ASSET_DIR.resolve("structure-viewer.js"));
        copyResource("/ui/vendor/3Dmol-min.js", ASSET_DIR.resolve("vendor/3Dmol-min.js"));
        return ASSET_DIR;
    }

    private static void copyResource(String resourcePath, Path target) throws IOException {
        try (InputStream input = ChromiumStructureViewer.class.getResourceAsStream(resourcePath)) {
            if (input == null) {
                throw new IOException("Missing viewer resource " + resourcePath);
            }
            Files.createDirectories(target.getParent());
            Files.copy(input, target, StandardCopyOption.REPLACE_EXISTING);
        }
    }

    private static String extractJsonValue(String json, String key) {
        String marker = "\"" + key + "\":\"";
        int start = json.indexOf(marker);
        if (start < 0) {
            return "";
        }
        int valueStart = start + marker.length();
        int valueEnd = json.indexOf('"', valueStart);
        if (valueEnd < 0) {
            return "";
        }
        return json.substring(valueStart, valueEnd)
                .replace("\\n", "\n")
                .replace("\\r", "\r")
                .replace("\\\"", "\"")
                .replace("\\\\", "\\");
    }

    private static String quoteJs(String value) {
        if (value == null) {
            return "null";
        }
        return "\"" + value
                .replace("\\", "\\\\")
                .replace("\"", "\\\"")
                .replace("\r", "\\r")
                .replace("\n", "\\n") + "\"";
    }

    private static String escapeHtml(String value) {
        if (value == null) {
            return "";
        }
        return value.replace("&", "&amp;")
                .replace("<", "&lt;")
                .replace(">", "&gt;");
    }
}
